import argparse
import json
import re
import sys
import urllib.error
import urllib.request
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

try:
    import pdfplumber
except ImportError:  # pragma: no cover
    pdfplumber = None

try:
    from pypdf import PdfReader
except ImportError:  # pragma: no cover
    PdfReader = None

try:
    import camelot
except ImportError:  # pragma: no cover
    camelot = None

try:
    import pandas as pd
except ImportError:  # pragma: no cover
    pd = None


COLUMNS = [
    "Country of Origin",
    "Title",
    "Authors",
    "Year",
    "Journal",
    "DOI",
    "Keywords",
    "Primary Organism",
    "Secondary Organism",
    "Gene Family",
    "Genome Source",
    "HMMER (with Pfam domains)",
    "BLAST (E-value cut-off used)",
    "SMART",
    "Other (Tools for Domain Analysis)",
    "CDD",
    "MEGA",
    "IQ-TREE",
    "MEME",
    "Others (Tools for Motif Analysis)",
    "MCScanX",
    "Other (Tools for Synteny)",
    "TBTools",
    "Other (Tools for Visualization)",
    "Experimental Validation (qPCR)",
    "Other Tools Used In Methodology",
    "Number of Genes Identified",
    "Phylogenetic Structure",
    "Gene Duplication Events (TANDEM/SEGMENTAL)",
    "Expression Analysis (RNA-seq)",
]

COUNTRIES = [
    "Pakistan", "India", "China", "Japan", "South Korea", "Bangladesh", "Thailand", "Vietnam", "Malaysia",
    "Indonesia", "Philippines", "Saudi Arabia", "United Arab Emirates", "Iran", "Turkey", "Egypt",
    "South Africa", "Nigeria", "Kenya", "United States", "USA", "Canada", "Mexico", "Brazil", "Argentina",
    "United Kingdom", "Germany", "France", "Italy", "Spain", "Netherlands", "Belgium", "Switzerland",
    "Sweden", "Norway", "Denmark", "Finland", "Poland", "Australia", "New Zealand",
]

JOURNAL_HINTS = [
    "Biochemical and Biophysical Research Communications",
    "Plant Gene",
    "Horticulturae",
    "Genes",
    "Genetics Research",
    "Journal of Biomolecular Structure and Dynamics",
    "Genomics & Informatics",
    "Genomics Inform",
    "BMC Genomics",
    "BMC Plant Biology",
    "Frontiers in Plant Science",
    "International Journal of Molecular Sciences",
    "Plant Physiology",
    "Plant Journal",
    "Scientific Reports",
    "Nucleic Acids Research",
]

COMMON_SPECIES = [
    "Arabidopsis thaliana", "Oryza sativa", "Zea mays", "Glycine max", "Solanum lycopersicum",
    "Capsicum annuum", "Triticum aestivum", "Hordeum vulgare", "Brassica napus", "Brassica rapa",
    "Gossypium hirsutum", "Malus domestica", "Populus trichocarpa", "Vitis vinifera",
    "Nicotiana tabacum", "Medicago truncatula", "Sorghum bicolor", "Setaria italica",
    "Camellia nitidissima", "Solanum tuberosum", "Bubalus bubalis",
]

KNOWN_GENE_FAMILIES = [
    "WRKY", "MYB", "NAC", "ERF", "AP2/ERF", "bHLH", "MADS", "HD-ZIP", "GRAS", "ARF", "TCP", "SBP",
    "GATA", "DOF", "HSF", "C2H2", "WOX", "SPL", "BBX", "PAL", "LIM", "Aquaporin", "FGF", "Trihelix",
    "TPL/TPR", "Lectin receptor-like kinase", "LecRLK", "MscS-like",
]

DOMAIN_TOOLS = ["HMMER", "Pfam", "SMART", "CDD", "InterPro", "InterProScan", "PROSITE", "ScanProsite"]
MOTIF_TOOLS = ["MEME", "FIMO", "Tomtom", "HOMER", "DREME", "STREME"]
SYNTENY_TOOLS = ["MCScanX", "JCVI", "SynMap", "CoGe", "WGDI", "i-ADHoRe"]
VIS_TOOLS = ["TBtools", "GSDS", "EvolView", "iTOL", "Circos", "MapChart"]
METHOD_TOOLS = [
    "BLAST", "BLASTP", "HMMER", "MAFFT", "MUSCLE", "ClustalW", "MEGA", "IQ-TREE", "MEME", "MCScanX",
    "TBtools", "PlantCARE", "TargetP", "WoLF PSORT", "InterProScan", "SMART", "CDD", "ExPASy",
]

SECTION_HEADERS = {
    "abstract": ["abstract"],
    "keywords": ["keywords", "key words"],
    "methods": ["materials and methods", "methods", "methodology", "experimental procedures"],
    "results": ["results"],
    "discussion": ["discussion"],
    "conclusion": ["conclusion", "conclusions"],
}

NOISE_TOKENS = [
    "license", "licensee", "copyright", "citation", "published:", "received", "accepted", "revised",
    "academic editor", "article history", "journal homepage", "contents lists available",
]

NUMBER_WORDS = {
    "one": 1, "two": 2, "three": 3, "four": 4, "five": 5, "six": 6, "seven": 7, "eight": 8, "nine": 9,
    "ten": 10, "eleven": 11, "twelve": 12, "thirteen": 13, "fourteen": 14, "fifteen": 15, "sixteen": 16,
    "seventeen": 17, "eighteen": 18, "nineteen": 19, "twenty": 20,
}


def ensure_dependencies() -> None:
    if pdfplumber is None and PdfReader is None:
        raise RuntimeError("Install `pdfplumber` or `pypdf` in your local environment.")
    if pd is None:
        raise RuntimeError("Install `pandas` in your local environment.")


def normalize_text(text: str) -> str:
    replacements = {
        "\u00a0": " ",
        "\ufb00": "ff",
        "\ufb01": "fi",
        "\ufb02": "fl",
        "\ufb03": "ffi",
        "\ufb04": "ffl",
        "\u2010": "-",
        "\u2011": "-",
        "\u2012": "-",
        "\u2013": "-",
        "\u2014": "-",
        "\u2212": "-",
        "\u00d7": "x",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    text = re.sub(r"(?<=[a-z])(?=[A-Z])", " ", text)
    text = re.sub(r"[ \t]+", " ", text)
    text = re.sub(r"\n{3,}", "\n\n", text)
    return text.strip()


def flatten(text: str) -> str:
    return re.sub(r"\s+", " ", normalize_text(text))


def clean_value(value: object) -> str:
    if value is None:
        return ""
    text = re.sub(r"\s+", " ", str(value)).strip()
    if text.lower() in {"nan", "none", "null"}:
        return ""
    return text


def first(items: Iterable[str]) -> str:
    for item in items:
        item = clean_value(item)
        if item:
            return item
    return ""


def dedupe(items: Iterable[str]) -> list[str]:
    seen = set()
    output = []
    for item in items:
        item = clean_value(item).strip(" ,;:.")
        if not item:
            continue
        key = item.lower()
        if key in seen:
            continue
        seen.add(key)
        output.append(item)
    return output


def numeric_to_int(token: str) -> Optional[int]:
    token = clean_value(token).lower().replace("-", " ")
    if token.isdigit():
        return int(token)
    parts = [part for part in token.split() if part]
    if not parts:
        return None
    total = 0
    for part in parts:
        value = NUMBER_WORDS.get(part)
        if value is None:
            return None
        total += value
    return total if total else None


@dataclass
class FieldResult:
    value: str
    confidence: float
    evidence: str
    source: str


@dataclass
class DocumentContext:
    pdf_path: Path
    full_text: str
    first_page_text: str
    table_text: str
    sections: dict[str, str]

    @property
    def combined_text(self) -> str:
        return normalize_text(self.full_text + ("\n\n" + self.table_text if self.table_text else ""))

    @property
    def methods_text(self) -> str:
        return "\n".join(
            part for part in [
                self.sections.get("methods", ""),
                self.sections.get("results", ""),
                self.sections.get("discussion", ""),
            ] if part
        ).strip()


class OllamaHelper:
    def __init__(self, model: str):
        self.model = model

    def available(self) -> bool:
        try:
            req = urllib.request.Request("http://127.0.0.1:11434/api/tags", method="GET")
            with urllib.request.urlopen(req, timeout=2) as response:
                payload = json.loads(response.read().decode("utf-8"))
            return any(model.get("name", "").startswith(self.model) for model in payload.get("models", []))
        except Exception:
            return False

    def ask(self, prompt: str) -> str:
        body = json.dumps(
            {"model": self.model, "prompt": prompt, "stream": False, "options": {"temperature": 0}},
        ).encode("utf-8")
        req = urllib.request.Request(
            "http://127.0.0.1:11434/api/generate",
            data=body,
            headers={"Content-Type": "application/json"},
            method="POST",
        )
        with urllib.request.urlopen(req, timeout=120) as response:
            payload = json.loads(response.read().decode("utf-8"))
        return clean_value(payload.get("response", ""))


class CompleteLocalExtractor:
    def __init__(
        self,
        review_threshold: float = 0.65,
        use_tables: bool = False,
        ollama_model: Optional[str] = None,
    ):
        self.review_threshold = review_threshold
        self.use_tables = use_tables
        self.ollama = OllamaHelper(ollama_model) if ollama_model else None
        self.ollama_enabled = bool(self.ollama and self.ollama.available())

    def reader_score(self, text: str) -> float:
        lines = [line.strip() for line in text.splitlines() if line.strip()]
        if not lines:
            return float("-inf")
        score = sum(len(re.findall(r"[A-Za-z]{2,}", line)) for line in lines[:80])
        score -= sum(1 for line in lines[:80] if len(line) <= 2) * 5
        score -= sum(1 for line in lines[:80] if re.fullmatch(r"[A-Za-z]", line)) * 8
        score += 12 if any("gene family" in line.lower() for line in lines[:20]) else 0
        score += 8 if any("abstract" == line.lower() for line in lines[:30]) else 0
        return score

    def extract_reader_outputs(self, pdf_path: Path) -> dict[str, str]:
        outputs: dict[str, str] = {}
        if pdfplumber is not None:
            try:
                with pdfplumber.open(str(pdf_path)) as pdf:
                    parts = [(page.extract_text() or "").strip() for page in pdf.pages]
                text = normalize_text("\n\n".join(part for part in parts if part))
                if text:
                    outputs["pdfplumber"] = text
            except Exception:
                pass
        if PdfReader is not None:
            try:
                reader = PdfReader(str(pdf_path))
                parts = [(page.extract_text() or "").strip() for page in reader.pages]
                text = normalize_text("\n\n".join(part for part in parts if part))
                if text:
                    outputs["pypdf"] = text
            except Exception:
                pass
        return outputs

    def extract_best_text(self, pdf_path: Path) -> str:
        outputs = self.extract_reader_outputs(pdf_path)
        if not outputs:
            return ""
        return max(outputs.values(), key=self.reader_score)

    def extract_best_first_page(self, pdf_path: Path) -> str:
        candidates = []
        if pdfplumber is not None:
            try:
                with pdfplumber.open(str(pdf_path)) as pdf:
                    if pdf.pages:
                        text = normalize_text(pdf.pages[0].extract_text() or "")
                        if text:
                            candidates.append(text)
            except Exception:
                pass
        if PdfReader is not None:
            try:
                reader = PdfReader(str(pdf_path))
                if reader.pages:
                    text = normalize_text(reader.pages[0].extract_text() or "")
                    if text:
                        candidates.append(text)
            except Exception:
                pass
        if not candidates:
            return ""
        return max(candidates, key=self.reader_score)

    def extract_table_text(self, pdf_path: Path) -> str:
        if not self.use_tables or camelot is None:
            return ""
        snippets = []
        for flavor in ("lattice", "stream"):
            try:
                tables = camelot.read_pdf(str(pdf_path), pages="all", flavor=flavor)
            except Exception:
                continue
            for table in tables:
                try:
                    for row in table.df.fillna("").values.tolist():
                        joined = " | ".join(cell.strip() for cell in row if str(cell).strip())
                        if joined:
                            snippets.append(joined)
                except Exception:
                    continue
            if snippets:
                break
        return normalize_text("\n".join(snippets))

    def lines(self, text: str, limit: Optional[int] = None) -> list[str]:
        result = [line.strip() for line in normalize_text(text).splitlines() if line.strip()]
        return result[:limit] if limit else result

    def sections(self, text: str) -> dict[str, str]:
        found = {name: [] for name in SECTION_HEADERS}
        current = None
        for line in self.lines(text):
            low = line.lower().strip(" :")
            switched = False
            for section_name, headers in SECTION_HEADERS.items():
                if low in headers:
                    current = section_name
                    switched = True
                    break
            if switched:
                continue
            if current:
                found[current].append(line)
        return {name: "\n".join(value).strip() for name, value in found.items()}

    def build_context(self, pdf_path: Path) -> DocumentContext:
        full_text = self.extract_best_text(pdf_path)
        first_page = self.extract_best_first_page(pdf_path) or full_text[:8000]
        table_text = self.extract_table_text(pdf_path)
        combined = normalize_text(full_text + ("\n\n" + table_text if table_text else ""))
        return DocumentContext(
            pdf_path=pdf_path,
            full_text=full_text,
            first_page_text=first_page,
            table_text=table_text,
            sections=self.sections(combined),
        )

    def evidence(self, text: str, patterns: Iterable[str]) -> str:
        for line in self.lines(text, 200):
            if any(pattern.lower() in line.lower() for pattern in patterns):
                return clean_value(line)[:220]
        return clean_value(first(self.lines(text, 20)))[:220]

    def looks_like_journal(self, line: str) -> bool:
        low = line.lower()
        return (
            line in JOURNAL_HINTS
            or "journal homepage" in low
            or re.search(r"\b(volume|issn|eissn|research article|original article)\b", low) is not None
            or any(hint.lower() in low for hint in ["plant gene", "genetics research", "genomics inform"])
        )

    def looks_like_author_line(self, line: str) -> bool:
        low = line.lower()
        if any(token in low for token in NOISE_TOKENS):
            return False
        if any(token in low for token in ["university", "department", "college", "institute", "faculty"]):
            return False
        if "@" in line or self.looks_like_journal(line):
            return False
        name_hits = len(re.findall(r"[A-Z][a-z]+(?:\s+[A-Z][a-z]+)+", line))
        return name_hits >= 2 or (";" in line and any(char.isalpha() for char in line))

    def clean_author_line(self, line: str) -> str:
        line = re.sub(r"(?<=[A-Za-z])\d+", "", line)
        line = re.sub(r"[*#\u2020†]+", "", line)
        line = re.sub(r"\s+", " ", line)
        return line.strip(" ,;")

    def extract_title(self, first_page_text: str) -> str:
        candidates = []
        front = self.lines(first_page_text, 35)
        for idx, line in enumerate(front):
            low = line.lower()
            if len(line) < 20 or len(line) > 240:
                continue
            if any(token in low for token in NOISE_TOKENS):
                continue
            if self.looks_like_author_line(line) or self.looks_like_journal(line):
                continue
            score = 0
            if idx <= 8:
                score += 5
            if re.search(r"\b(genome-wide|gene family|identification|analysis|characterization|evolutionary)\b", low):
                score += 6
            if line[:1].isupper():
                score += 2
            if re.search(r"[.!?]$", line):
                score -= 2
            if re.search(r"\b(volume|issue|issn|eissn)\b", low):
                score -= 4
            candidates.append((score, line, idx))
        if not candidates:
            return ""
        candidates.sort(key=lambda item: item[0], reverse=True)
        best_line = candidates[0][1]
        best_idx = candidates[0][2]
        merged = [best_line]
        for line in front[best_idx + 1: best_idx + 4]:
            if self.looks_like_author_line(line) or self.looks_like_journal(line):
                break
            if any(token in line.lower() for token in NOISE_TOKENS):
                break
            if len(line) > 120:
                break
            if line[:1].islower() or re.search(r"\([A-Z][a-z]+\s+[a-z]+", line):
                merged.append(line)
            else:
                break
        text = " ".join(merged)
        text = re.sub(r"\bin([A-Z])", r"in \1", text)
        text = re.sub(r"\bof([A-Z])", r"of \1", text)
        return clean_value(text)

    def extract_authors(self, first_page_text: str, title: str) -> str:
        front = self.lines(first_page_text, 40)
        for idx, line in enumerate(front):
            if title and clean_value(line) == clean_value(title):
                for follow in front[idx + 1: idx + 8]:
                    if self.looks_like_author_line(follow):
                        return self.clean_author_line(follow)
            if self.looks_like_author_line(line):
                return self.clean_author_line(line)
        return ""

    def extract_year(self, text: str) -> str:
        years = re.findall(r"\b(19\d{2}|20\d{2})\b", text[:10000])
        valid = [year for year in years if 1990 <= int(year) <= 2035]
        return Counter(valid).most_common(1)[0][0] if valid else ""

    def extract_journal(self, first_page_text: str) -> str:
        front = self.lines(first_page_text, 25)
        for line in front:
            if line in JOURNAL_HINTS:
                return line
        for line in front:
            match = re.search(r"\.([A-Z][A-Za-z& \-]+)\s+(19|20)\d{2}\b", line)
            if match:
                return clean_value(match.group(1))
        for line in front:
            if self.looks_like_journal(line) and line.lower() not in {"article", "research article", "original article"}:
                return clean_value(line)
        return ""

    def extract_doi(self, text: str) -> str:
        match = re.search(r"\b10\.\d{4,9}/[-._;()/:A-Z0-9]+\b", text, flags=re.IGNORECASE)
        return clean_value(match.group(0)) if match else ""

    def extract_keywords(self, ctx: DocumentContext) -> str:
        keywords = ctx.sections.get("keywords", "")
        if not keywords:
            match = re.search(r"keywords?\s*[:\-]\s*(.{5,250})", ctx.combined_text, flags=re.IGNORECASE)
            keywords = match.group(1) if match else ""
        keywords = re.split(r"\n|\. ", keywords)[0]
        if "identified" in keywords.lower() or len(re.findall(r"[A-Za-z]{4,}", keywords)) > 18:
            return ""
        return ", ".join(dedupe(re.split(r"[;,]", keywords)[:15]))

    def extract_country(self, first_page_text: str) -> str:
        found = []
        for country in COUNTRIES:
            if re.search(r"\b" + re.escape(country) + r"\b", first_page_text, flags=re.IGNORECASE):
                found.append("United States" if country == "USA" else country)
        return Counter(found).most_common(1)[0][0] if found else ""

    def extract_organisms(self, title: str, abstract: str, keywords: str) -> tuple[str, str]:
        context = "\n".join([title, abstract, keywords])
        found = []
        for species in COMMON_SPECIES:
            if species.lower() in context.lower():
                found.extend([species] * 3)
        found.extend(re.findall(r"\b([A-Z][a-z]+\s+[a-z][a-z\-]+)\b", context))
        filtered = []
        for item in found:
            if len(item.split()) != 2:
                continue
            low = item.lower()
            if low.startswith(("this study", "many studies", "overexpression of")):
                continue
            filtered.append(item)
        common = [name for name, _ in Counter(filtered).most_common(2)]
        return (common[0] if common else "", common[1] if len(common) > 1 else "")

    def extract_gene_family(self, title: str, abstract: str) -> str:
        for sample in [title, abstract, title + "\n" + abstract]:
            for pattern in [
                r"\b([A-Za-z0-9/\-]+)\s+gene family\b",
                r"\b([A-Za-z0-9/\-]+)\s+transcription factor(?:s)?\b",
                r"\b([A-Za-z0-9/\-]+)\s+family genes\b",
                r"\b([A-Za-z0-9/\-]+)\s+like gene family\b",
            ]:
                match = re.search(pattern, sample, flags=re.IGNORECASE)
                if match:
                    value = re.sub(r"^(the|a|an)", "", match.group(1), flags=re.IGNORECASE).strip()
                    return clean_value(value)
        for family in KNOWN_GENE_FAMILIES:
            if re.search(r"\b" + re.escape(family) + r"\b", title + "\n" + abstract, flags=re.IGNORECASE):
                return family
        return ""

    def extract_genome_source(self, text: str) -> str:
        sources = ["Phytozome", "NCBI", "Ensembl", "EnsemblPlants", "TAIR", "Gramene", "JGI", "MaizeGDB", "SGN"]
        return ", ".join(dedupe(source for source in sources if source.lower() in text.lower()))

    def extract_hmmer(self, text: str) -> str:
        used = re.search(r"\b(HMMER|hmmsearch|hmmscan)\b", text, flags=re.IGNORECASE) is not None
        pfams = dedupe(re.findall(r"\bPF\d{5}\b", text))
        has_pfam = re.search(r"\bPfam\b", text, flags=re.IGNORECASE) is not None
        if not used and not pfams and not has_pfam:
            return "No"
        if pfams:
            return "Yes; Pfam domains: " + ", ".join(pfams[:15])
        if has_pfam:
            return "Yes; Pfam mentioned"
        return "Yes"

    def extract_blast(self, text: str) -> str:
        if re.search(r"\bBLAST[PNX]?\b", text, flags=re.IGNORECASE) is None:
            return "No"
        flat = flatten(text)
        for pattern in [
            r"BLAST[PNX]?\b.{0,80}?E-?value(?:\s+cut-?off|\s+threshold)?(?:\s+of)?\s*[:=<>~]?\s*([0-9.]+(?:e-?\d+)?)",
            r"E-?value(?:\s+cut-?off|\s+threshold)?(?:\s+of)?\s*[:=<>~]?\s*([0-9.]+(?:e-?\d+)?)",
            r"([0-9.]+\s*x\s*10-?\d+)",
        ]:
            match = re.search(pattern, flat, flags=re.IGNORECASE)
            if match:
                value = match.group(1).replace(" ", "").replace("x10", "e").replace("10-", "e-")
                return value
        return "Yes"

    def tool_yes_no(self, text: str, aliases: Iterable[str]) -> str:
        return "Yes" if any(re.search(r"\b" + re.escape(alias) + r"\b", text, flags=re.IGNORECASE) for alias in aliases) else "No"

    def extract_other_tools(self, text: str, pool: Iterable[str], exclude: Iterable[str]) -> str:
        exclude_set = {item.lower() for item in exclude}
        return ", ".join(dedupe(tool for tool in pool if tool.lower() not in exclude_set and tool.lower() in text.lower()))

    def extract_qpcr(self, text: str) -> str:
        return "Yes" if re.search(r"\b(qPCR|RT-qPCR|qRT-PCR|real-time PCR|SYBR Green)\b", text, flags=re.IGNORECASE) else "No"

    def extract_gene_count(self, title: str, abstract: str, results: str) -> str:
        search_space = flatten("\n".join([title, abstract, results]))
        patterns = [
            r"\ba total of ((?:\d{1,4})|(?:[A-Za-z\-]+)) .*? genes? (?:were )?identified\b",
            r"\bidentified ((?:\d{1,4})|(?:[A-Za-z\-]+)) .*? genes?\b",
            r"\bfound ((?:\d{1,4})|(?:[A-Za-z\-]+)) .*? genes?\b",
            r"\b((?:\d{1,4})|(?:[A-Za-z\-]+)) members of the .*? gene family\b",
        ]
        for pattern in patterns:
            match = re.search(pattern, search_space, flags=re.IGNORECASE)
            if match:
                raw = match.group(1)
                if raw.isdigit():
                    return raw
                converted = numeric_to_int(raw)
                if converted is not None:
                    return str(converted)
        return ""

    def extract_phylo_structure(self, text: str) -> str:
        for pattern in [
            r"\b(\d+)\s+(clades|groups|subfamilies|clusters)\b",
            r"\bclassified into (\d+)\s+(clades|groups|subfamilies|clusters)\b",
        ]:
            match = re.search(pattern, flatten(text), flags=re.IGNORECASE)
            if match:
                return f"{match.group(1)} {match.group(2)}"
        return ""

    def extract_duplication(self, text: str) -> str:
        events = []
        if re.search(r"\btandem duplication\b", text, flags=re.IGNORECASE):
            events.append("TANDEM")
        if re.search(r"\bsegmental duplication\b", text, flags=re.IGNORECASE):
            events.append("SEGMENTAL")
        return "/".join(events)

    def extract_rnaseq(self, text: str) -> str:
        return "Yes" if re.search(r"\b(RNA-seq|RNA sequencing|transcriptome|FPKM|TPM)\b", text, flags=re.IGNORECASE) else "No"

    def confidence(self, field: str, value: str) -> float:
        value = clean_value(value)
        if not value:
            return 0.0
        if field == "DOI":
            return 0.98 if re.match(r"^10\.\d{4,9}/", value, flags=re.IGNORECASE) else 0.3
        if field == "Year":
            return 0.95 if re.fullmatch(r"(19|20)\d{2}", value) else 0.25
        if field in {"SMART", "CDD", "MEGA", "IQ-TREE", "MEME", "MCScanX", "TBTools", "Experimental Validation (qPCR)", "Expression Analysis (RNA-seq)"}:
            return 0.9 if value in {"Yes", "No"} else 0.4
        if field == "Number of Genes Identified":
            return 0.9 if re.fullmatch(r"\d{1,4}", value) else 0.3
        score = 0.55
        if len(value) >= 5:
            score += 0.1
        if field in {"Title", "Authors", "Journal"}:
            score += 0.1
        if re.search(r"\b(revised|published|citation|license|editor)\b", value, flags=re.IGNORECASE):
            score -= 0.45
        if field in {"Primary Organism", "Secondary Organism"} and len(value.split()) != 2:
            score -= 0.2
        if field == "Journal" and re.search(r"\bissn\b", value, flags=re.IGNORECASE):
            score -= 0.3
        return round(max(0.0, min(1.0, score)), 2)

    def make_result(self, field: str, value: str, evidence_text: str, source: str, patterns: Iterable[str]) -> FieldResult:
        return FieldResult(
            value=clean_value(value),
            confidence=self.confidence(field, value),
            evidence=self.evidence(evidence_text, patterns),
            source=source,
        )

    def llm_repair(self, field: str, ctx: DocumentContext, current: FieldResult) -> FieldResult:
        if not self.ollama_enabled or current.confidence >= self.review_threshold:
            return current
        evidence_block = "\n\n".join(
            part for part in [
                ctx.first_page_text[:3000],
                ctx.sections.get("abstract", "")[:3000],
                ctx.methods_text[:3000],
            ] if part
        )
        prompt = (
            f"Extract only the value for this field from the paper text.\n"
            f"Field: {field}\n"
            f"Return one short value only. If missing, return EMPTY.\n\n"
            f"Paper text:\n{evidence_block}"
        )
        try:
            answer = self.ollama.ask(prompt)
        except (urllib.error.URLError, TimeoutError, OSError, ValueError):
            return current
        if not answer or answer.upper() == "EMPTY":
            return current
        repaired = self.make_result(field, answer, evidence_block, "ollama", [field])
        return repaired if repaired.confidence > current.confidence else current

    def extract_fields(self, ctx: DocumentContext) -> dict[str, FieldResult]:
        title = self.extract_title(ctx.first_page_text)
        authors = self.extract_authors(ctx.first_page_text, title)
        year = self.extract_year(ctx.first_page_text or ctx.full_text)
        journal = self.extract_journal(ctx.first_page_text)
        doi = self.extract_doi(ctx.combined_text)
        keywords = self.extract_keywords(ctx)
        abstract = ctx.sections.get("abstract", "")
        primary, secondary = self.extract_organisms(title, abstract, keywords)
        gene_family = self.extract_gene_family(title, abstract)
        genome_source = self.extract_genome_source(ctx.combined_text)
        methods_text = ctx.methods_text or ctx.combined_text

        data = {
            "Country of Origin": self.make_result("Country of Origin", self.extract_country(ctx.first_page_text), ctx.first_page_text, "front", COUNTRIES),
            "Title": self.make_result("Title", title, ctx.first_page_text, "front", [title]),
            "Authors": self.make_result("Authors", authors, ctx.first_page_text, "front", [authors]),
            "Year": self.make_result("Year", year, ctx.first_page_text, "front", ["received", "published", "accepted"]),
            "Journal": self.make_result("Journal", journal, ctx.first_page_text, "front", JOURNAL_HINTS),
            "DOI": self.make_result("DOI", doi, ctx.combined_text, "full", ["doi", "https://doi.org"]),
            "Keywords": self.make_result("Keywords", keywords, ctx.sections.get("keywords", "") or ctx.combined_text, "keywords", ["keywords"]),
            "Primary Organism": self.make_result("Primary Organism", primary, "\n".join([title, abstract, keywords]), "context", COMMON_SPECIES),
            "Secondary Organism": self.make_result("Secondary Organism", secondary, "\n".join([title, abstract, keywords]), "context", COMMON_SPECIES),
            "Gene Family": self.make_result("Gene Family", gene_family, "\n".join([title, abstract]), "front+abstract", KNOWN_GENE_FAMILIES),
            "Genome Source": self.make_result("Genome Source", genome_source, ctx.combined_text, "full", ["Phytozome", "NCBI", "Ensembl", "TAIR", "JGI"]),
            "HMMER (with Pfam domains)": self.make_result("HMMER (with Pfam domains)", self.extract_hmmer(ctx.combined_text), ctx.combined_text, "full", ["HMMER", "Pfam"]),
            "BLAST (E-value cut-off used)": self.make_result("BLAST (E-value cut-off used)", self.extract_blast(ctx.combined_text), ctx.combined_text, "full", ["BLAST", "E-value"]),
            "SMART": self.make_result("SMART", self.tool_yes_no(ctx.combined_text, ["SMART"]), ctx.combined_text, "full", ["SMART"]),
            "Other (Tools for Domain Analysis)": self.make_result("Other (Tools for Domain Analysis)", self.extract_other_tools(ctx.combined_text, DOMAIN_TOOLS, ["HMMER", "Pfam", "SMART", "CDD"]), ctx.combined_text, "full", DOMAIN_TOOLS),
            "CDD": self.make_result("CDD", self.tool_yes_no(ctx.combined_text, ["CDD", "Conserved Domain Database", "NCBI-CDD"]), ctx.combined_text, "full", ["CDD"]),
            "MEGA": self.make_result("MEGA", self.tool_yes_no(ctx.combined_text, ["MEGA"]), ctx.combined_text, "full", ["MEGA"]),
            "IQ-TREE": self.make_result("IQ-TREE", self.tool_yes_no(ctx.combined_text, ["IQ-TREE", "IQTREE", "iqtree"]), ctx.combined_text, "full", ["IQ-TREE"]),
            "MEME": self.make_result("MEME", self.tool_yes_no(ctx.combined_text, ["MEME"]), ctx.combined_text, "full", ["MEME"]),
            "Others (Tools for Motif Analysis)": self.make_result("Others (Tools for Motif Analysis)", self.extract_other_tools(ctx.combined_text, MOTIF_TOOLS, ["MEME"]), ctx.combined_text, "full", MOTIF_TOOLS),
            "MCScanX": self.make_result("MCScanX", self.tool_yes_no(ctx.combined_text, ["MCScanX"]), ctx.combined_text, "full", ["MCScanX"]),
            "Other (Tools for Synteny)": self.make_result("Other (Tools for Synteny)", self.extract_other_tools(ctx.combined_text, SYNTENY_TOOLS, ["MCScanX"]), ctx.combined_text, "full", SYNTENY_TOOLS),
            "TBTools": self.make_result("TBTools", self.tool_yes_no(ctx.combined_text, ["TBtools", "TBTools"]), ctx.combined_text, "full", ["TBtools"]),
            "Other (Tools for Visualization)": self.make_result("Other (Tools for Visualization)", self.extract_other_tools(ctx.combined_text, VIS_TOOLS, ["TBtools"]), ctx.combined_text, "full", VIS_TOOLS),
            "Experimental Validation (qPCR)": self.make_result("Experimental Validation (qPCR)", self.extract_qpcr(methods_text), methods_text, "methods", ["qPCR", "RT-qPCR"]),
            "Other Tools Used In Methodology": self.make_result("Other Tools Used In Methodology", self.extract_other_tools(methods_text, METHOD_TOOLS, []), methods_text, "methods", METHOD_TOOLS),
            "Number of Genes Identified": self.make_result("Number of Genes Identified", self.extract_gene_count(title, abstract, ctx.sections.get("results", "")), "\n".join([title, abstract, ctx.sections.get("results", "")]), "title+abstract+results", ["identified", "genes"]),
            "Phylogenetic Structure": self.make_result("Phylogenetic Structure", self.extract_phylo_structure(ctx.combined_text), ctx.combined_text, "full", ["clades", "groups", "subfamilies"]),
            "Gene Duplication Events (TANDEM/SEGMENTAL)": self.make_result("Gene Duplication Events (TANDEM/SEGMENTAL)", self.extract_duplication(ctx.combined_text), ctx.combined_text, "full", ["tandem duplication", "segmental duplication"]),
            "Expression Analysis (RNA-seq)": self.make_result("Expression Analysis (RNA-seq)", self.extract_rnaseq(methods_text), methods_text, "methods", ["RNA-seq", "transcriptome"]),
        }
        if self.ollama_enabled:
            for field in ["Title", "Authors", "Journal", "Primary Organism", "Secondary Organism", "Gene Family", "Number of Genes Identified"]:
                data[field] = self.llm_repair(field, ctx, data[field])
        return data

    def flatten_row(self, fields: dict[str, FieldResult]) -> dict[str, str]:
        return {column: fields[column].value for column in COLUMNS}

    def review_rows(self, pdf_name: str, fields: dict[str, FieldResult]) -> list[dict]:
        rows = []
        for field, result in fields.items():
            if result.confidence < self.review_threshold:
                rows.append(
                    {
                        "paper": pdf_name,
                        "field": field,
                        "value": result.value,
                        "confidence": result.confidence,
                        "source": result.source,
                        "evidence": result.evidence,
                    }
                )
        return rows

    def summary_row(self, pdf_name: str, fields: dict[str, FieldResult]) -> dict:
        confidences = [result.confidence for result in fields.values()]
        low = [field for field, result in fields.items() if result.confidence < self.review_threshold]
        return {
            "paper": pdf_name,
            "average_confidence": round(sum(confidences) / len(confidences), 2),
            "low_confidence_fields": ", ".join(low),
            "needs_review": "Yes" if low else "No",
        }

    def process_directory(self, input_dir: Path, output_dir: Path, run_name: str, excel: bool) -> None:
        pdf_files = sorted(input_dir.rglob("*.pdf"))
        if not pdf_files:
            raise FileNotFoundError(f"No PDF files found in {input_dir}")

        output_dir.mkdir(parents=True, exist_ok=True)
        main_csv = output_dir / f"{run_name}.csv"
        review_csv = output_dir / f"{run_name}_review_queue.csv"
        summary_csv = output_dir / f"{run_name}_paper_summary.csv"
        evidence_jsonl = output_dir / f"{run_name}_evidence.jsonl"
        readme = output_dir / f"{run_name}_README.txt"

        rows = []
        reviews = []
        summaries = []
        evidence_lines = []

        print(f"Input folder: {input_dir}")
        print(f"Output folder: {output_dir}")
        print(f"Ollama repair enabled: {'Yes' if self.ollama_enabled else 'No'}")
        print(f"Table extraction enabled: {'Yes' if self.use_tables else 'No'}")

        for idx, pdf_path in enumerate(pdf_files, start=1):
            print(f"[{idx}/{len(pdf_files)}] Processing {pdf_path.name}")
            try:
                ctx = self.build_context(pdf_path)
                fields = self.extract_fields(ctx)
                rows.append(self.flatten_row(fields))
                reviews.extend(self.review_rows(pdf_path.name, fields))
                summaries.append(self.summary_row(pdf_path.name, fields))
                evidence_lines.append(
                    json.dumps(
                        {
                            "paper": pdf_path.name,
                            "fields": {
                                field: {
                                    "value": result.value,
                                    "confidence": result.confidence,
                                    "evidence": result.evidence,
                                    "source": result.source,
                                }
                                for field, result in fields.items()
                            },
                        },
                        ensure_ascii=False,
                    )
                )
            except Exception as exc:
                print(f"  Failed: {exc}")
                rows.append({column: "" for column in COLUMNS})
                summaries.append(
                    {
                        "paper": pdf_path.name,
                        "average_confidence": 0.0,
                        "low_confidence_fields": "ALL_FIELDS",
                        "needs_review": "Yes",
                    }
                )

        df = pd.DataFrame(rows, columns=COLUMNS)
        df.to_csv(main_csv, index=False, encoding="utf-8-sig")
        pd.DataFrame(reviews, columns=["paper", "field", "value", "confidence", "source", "evidence"]).to_csv(
            review_csv, index=False, encoding="utf-8-sig"
        )
        pd.DataFrame(summaries, columns=["paper", "average_confidence", "low_confidence_fields", "needs_review"]).to_csv(
            summary_csv, index=False, encoding="utf-8-sig"
        )
        evidence_jsonl.write_text("\n".join(evidence_lines), encoding="utf-8")

        readme.write_text(
            "\n".join(
                [
                    "Complete Local Research Paper Extractor",
                    "======================================",
                    f"Input folder: {input_dir}",
                    f"Total PDFs: {len(pdf_files)}",
                    f"Main CSV: {main_csv.name}",
                    f"Review queue: {review_csv.name}",
                    f"Paper summary: {summary_csv.name}",
                    f"Evidence JSONL: {evidence_jsonl.name}",
                    f"Review threshold: {self.review_threshold}",
                    f"Table extraction enabled: {'Yes' if self.use_tables else 'No'}",
                    f"Ollama repair enabled: {'Yes' if self.ollama_enabled else 'No'}",
                    "",
                    "Recommended workflow:",
                    "1. Open the main CSV for first-pass extraction.",
                    "2. Open the paper summary to find the weakest papers.",
                    "3. Open the review queue to manually verify only low-confidence cells.",
                    "4. Use the evidence JSONL if you want traceable field-by-field evidence.",
                ]
            ),
            encoding="utf-8",
        )

        if excel:
            xlsx_path = output_dir / f"{run_name}.xlsx"
            df.to_excel(xlsx_path, index=False)
            print(f"Saved Excel: {xlsx_path}")

        print(f"Saved CSV: {main_csv}")
        print(f"Saved review queue: {review_csv}")
        print(f"Saved paper summary: {summary_csv}")
        print(f"Saved evidence JSONL: {evidence_jsonl}")
        print(f"Saved run guide: {readme}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Free local research-paper extractor for genomics PDFs.")
    parser.add_argument("--input", default=".", help="Folder containing PDF files.")
    parser.add_argument("--output-dir", default="extractor_output", help="Folder where results will be written.")
    parser.add_argument("--run-name", default="local_extraction", help="Base name used for all output files.")
    parser.add_argument("--excel", action="store_true", help="Also create an XLSX copy of the main CSV.")
    parser.add_argument("--review-threshold", type=float, default=0.65, help="Low-confidence cutoff. Default: 0.65")
    parser.add_argument("--use-tables", action="store_true", help="Enable Camelot table extraction. Slower and sometimes noisier.")
    parser.add_argument("--ollama-model", default="", help="Optional local Ollama model name for low-confidence repair, e.g. qwen2.5:7b")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    ensure_dependencies()
    extractor = CompleteLocalExtractor(
        review_threshold=args.review_threshold,
        use_tables=args.use_tables,
        ollama_model=clean_value(args.ollama_model) or None,
    )
    extractor.process_directory(
        input_dir=Path(args.input).expanduser().resolve(),
        output_dir=Path(args.output_dir).expanduser().resolve(),
        run_name=clean_value(args.run_name) or "local_extraction",
        excel=args.excel,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
