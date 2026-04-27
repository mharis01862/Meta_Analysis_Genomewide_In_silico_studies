"""
genomewide_pdf_extractor.py
=====================
Genome-Wide Gene Family Study PDF Extractor
A fully local, API-free tool for extracting structured metadata from
genome-wide gene family identification study PDFs.

MIT License
Copyright (c) 2025

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Usage
-----
    python genomewide_pdf_extractor.py --input ./pdfs --output-dir ./results

Dependencies (install via pip)
------------------------------
    pip install pdfplumber pypdf pandas openpyxl

Optional (for Camelot table extraction – slower):
    pip install camelot-py[cv] ghostscript

Optional (for local LLM repair via Ollama):
    Install Ollama from https://ollama.com and pull a model, e.g.:
    ollama pull qwen2.5:7b
    Then pass: --ollama-model qwen2.5:7b
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import urllib.error
import urllib.request
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

# ---------------------------------------------------------------------------
# Optional dependency imports – graceful fallback if not installed
# ---------------------------------------------------------------------------

try:
    import pdfplumber
except ImportError:
    pdfplumber = None  # type: ignore[assignment]

try:
    from pypdf import PdfReader
except ImportError:
    PdfReader = None  # type: ignore[assignment]

try:
    import camelot
except ImportError:
    camelot = None  # type: ignore[assignment]

try:
    import pandas as pd
except ImportError:
    pd = None  # type: ignore[assignment]


# ===========================================================================
# DOMAIN KNOWLEDGE: column schema, controlled vocabularies, known entities
# ===========================================================================

#: All output columns, in order.
COLUMNS: list[str] = [
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

#: Recognised country names (used for affiliation extraction).
COUNTRIES: list[str] = [
    "Pakistan", "India", "China", "Japan", "South Korea", "Bangladesh",
    "Thailand", "Vietnam", "Malaysia", "Indonesia", "Philippines",
    "Saudi Arabia", "United Arab Emirates", "Iran", "Turkey", "Egypt",
    "South Africa", "Nigeria", "Kenya", "United States", "USA",
    "Canada", "Mexico", "Brazil", "Argentina", "United Kingdom",
    "Germany", "France", "Italy", "Spain", "Netherlands", "Belgium",
    "Switzerland", "Sweden", "Norway", "Denmark", "Finland", "Poland",
    "Australia", "New Zealand", "Russia", "Ukraine", "Ethiopia",
    "Tanzania", "Morocco", "Algeria", "Tunisia", "Ghana",
]

#: Journal names/hints for recognition on the first page.
JOURNAL_HINTS: list[str] = [
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
    "The Plant Journal",
    "Scientific Reports",
    "Nucleic Acids Research",
    "Plant and Cell Physiology",
    "New Phytologist",
    "Molecular Plant",
    "Planta",
    "Journal of Experimental Botany",
    "Plant Science",
    "Plant Molecular Biology",
    "Plant Cell Reports",
    "Plants",
    "Agronomy",
    "Agriculture",
    "Crop Science",
    "Gene",
    "Genomics",
    "Plant Cell",
    "Plant Physiology and Biochemistry",
    "Theoretical and Applied Genetics",
    "Tree Genetics & Genomes",
    "Forest Genetics and Tree Breeding",
]

#: Common plant/animal species studied in genome-wide gene family papers.
COMMON_SPECIES: list[str] = [
    "Arabidopsis thaliana", "Oryza sativa", "Zea mays", "Glycine max",
    "Solanum lycopersicum", "Capsicum annuum", "Triticum aestivum",
    "Hordeum vulgare", "Brassica napus", "Brassica rapa",
    "Gossypium hirsutum", "Gossypium arboreum", "Gossypium barbadense",
    "Malus domestica", "Populus trichocarpa", "Vitis vinifera",
    "Nicotiana tabacum", "Medicago truncatula", "Sorghum bicolor",
    "Setaria italica", "Camellia sinensis", "Camellia nitidissima",
    "Solanum tuberosum", "Bubalus bubalis", "Beta vulgaris",
    "Musa acuminata", "Cucumis sativus", "Cucurbita pepo",
    "Citrullus lanatus", "Phaseolus vulgaris", "Vigna radiata",
    "Vigna unguiculata", "Lotus japonicus", "Cajanus cajan",
    "Arachis hypogaea", "Helianthus annuus", "Sesamum indicum",
    "Linum usitatissimum", "Manihot esculenta", "Jatropha curcas",
    "Populus alba", "Eucalyptus grandis", "Pinus taeda",
    "Picea abies", "Phoenix dactylifera", "Elaeis guineensis",
    "Ananas comosus", "Saccharum officinarum",
    "Chenopodium quinoa", "Amaranthus hypochondriacus",
    "Spinacia oleracea", "Lactuca sativa",
    "Solanum melongena", "Ipomoea batatas",
    "Daucus carota", "Coriandrum sativum",
    "Panicum virgatum", "Brachypodium distachyon",
]

#: Well-known transcription factor / gene families.
KNOWN_GENE_FAMILIES: list[str] = [
    "WRKY", "MYB", "NAC", "ERF", "AP2/ERF", "AP2", "bHLH", "MADS",
    "HD-ZIP", "GRAS", "ARF", "TCP", "SBP", "GATA", "DOF", "HSF",
    "C2H2", "WOX", "SPL", "BBX", "PAL", "LIM", "Aquaporin", "FGF",
    "Trihelix", "TPL/TPR", "LecRLK", "MscS-like", "NBS-LRR",
    "LRR-RLK", "WAK", "GRF", "GIF", "TALE", "ZF-HD", "LBD",
    "EIN3/EIL", "ABI3/VP1", "bZIP", "RAV", "YABBY", "LATERAL ORGAN",
    "LOB", "ARR", "DREB", "CBF", "AREB", "TGA", "MYC", "JAZ",
    "BZR", "BES1", "HLH", "Dof", "R2R3-MYB", "R3-MYB",
    "Type-A ARR", "Type-B ARR", "RR", "CRF", "B3",
    "Ethylene Response Factor", "Dehydrin", "LEA",
    "Lipid Transfer Protein", "LTP", "Glutathione S-transferase",
    "GST", "Peroxidase", "Laccase", "Expansin", "Cellulose Synthase",
    "CesA", "Dirigent", "CYP450", "P450", "ABC transporter",
    "MATE", "NPF", "NRT", "PHT", "SUT", "SWEET", "PIN", "AUX",
    "LAX", "NHX", "HKT", "KT", "HAK", "AKT",
    "Receptor-Like Kinase", "RLK", "MAPK", "CDPK",
    "14-3-3", "HSP", "DnaJ", "Ubiquitin", "RING",
    "F-box", "SKP", "Cullin", "BTB", "PP2C",
]

#: Tool lists grouped by analytical function.
DOMAIN_TOOLS: list[str] = [
    "HMMER", "Pfam", "SMART", "CDD", "InterPro", "InterProScan",
    "PROSITE", "ScanProsite", "Superfamily", "PANTHER",
]
MOTIF_TOOLS: list[str] = [
    "MEME", "FIMO", "Tomtom", "HOMER", "DREME", "STREME",
    "MAST", "MCAST",
]
SYNTENY_TOOLS: list[str] = [
    "MCScanX", "JCVI", "SynMap", "CoGe", "WGDI", "i-ADHoRe",
    "MCScan", "SyMAP", "SynFind",
]
VIS_TOOLS: list[str] = [
    "TBtools", "GSDS", "EvolView", "iTOL", "Circos", "MapChart",
    "Cytoscape", "BioRender", "Inkscape", "R ggplot", "MEGA",
    "PlantGenIE",
]
METHOD_TOOLS: list[str] = [
    "BLAST", "BLASTP", "HMMER", "MAFFT", "MUSCLE", "ClustalW",
    "MEGA", "IQ-TREE", "MEME", "MCScanX", "TBtools", "PlantCARE",
    "TargetP", "WoLF PSORT", "InterProScan", "SMART", "CDD",
    "ExPASy", "SignalP", "TMHMM", "MODELLER", "PyMOL", "Swiss-Model",
    "PAML", "RAxML", "FastTree", "MISA", "Primer3",
]

#: Section-header aliases for structure detection.
SECTION_HEADERS: dict[str, list[str]] = {
    "abstract": ["abstract"],
    "keywords": ["keywords", "key words"],
    "introduction": ["introduction", "background"],
    "methods": [
        "materials and methods", "methods", "methodology",
        "experimental procedures", "materials & methods",
        "methods and materials",
    ],
    "results": ["results", "results and discussion"],
    "discussion": ["discussion"],
    "conclusion": ["conclusion", "conclusions", "concluding remarks"],
}

#: Substrings that indicate a line is metadata noise (not content).
NOISE_TOKENS: list[str] = [
    "license", "licensee", "copyright", "citation", "published:",
    "received", "accepted", "revised", "academic editor",
    "article history", "journal homepage", "contents lists available",
    "open access", "creative commons", "elsevier", "springer",
    "wiley", "mdpi", "frontiers in", "preprint",
]

#: Word → integer mapping for textual numbers in gene counts.
NUMBER_WORDS: dict[str, int] = {
    "one": 1, "two": 2, "three": 3, "four": 4, "five": 5,
    "six": 6, "seven": 7, "eight": 8, "nine": 9, "ten": 10,
    "eleven": 11, "twelve": 12, "thirteen": 13, "fourteen": 14,
    "fifteen": 15, "sixteen": 16, "seventeen": 17, "eighteen": 18,
    "nineteen": 19, "twenty": 20, "thirty": 30, "forty": 40,
    "fifty": 50, "sixty": 60, "seventy": 70, "eighty": 80,
    "ninety": 90, "hundred": 100,
}


# ===========================================================================
# UTILITY HELPERS
# ===========================================================================

def ensure_dependencies() -> None:
    """Raise helpful RuntimeError if critical libraries are missing."""
    if pdfplumber is None and PdfReader is None:
        raise RuntimeError(
            "Install at least one PDF reader:\n"
            "  pip install pdfplumber\n"
            "  pip install pypdf"
        )
    if pd is None:
        raise RuntimeError("Install pandas:\n  pip install pandas openpyxl")


def normalize_text(text: str) -> str:
    """Fix common PDF encoding artifacts and excessive whitespace."""
    _LIGATURES = {
        "\u00a0": " ",   # non-breaking space
        "\ufb00": "ff",  "\ufb01": "fi",  "\ufb02": "fl",
        "\ufb03": "ffi", "\ufb04": "ffl",
        "\u2010": "-",   "\u2011": "-",   "\u2012": "-",
        "\u2013": "-",   "\u2014": "-",   "\u2212": "-",
        "\u00d7": "x",   "\u2019": "'",   "\u2018": "'",
        "\u201c": '"',   "\u201d": '"',
    }
    for old, new in _LIGATURES.items():
        text = text.replace(old, new)
    # Merge camelCase splits introduced by some PDF renderers
    text = re.sub(r"(?<=[a-z])(?=[A-Z])", " ", text)
    text = re.sub(r"[ \t]+", " ", text)
    text = re.sub(r"\n{3,}", "\n\n", text)
    return text.strip()


def flatten(text: str) -> str:
    """Collapse all whitespace to a single space."""
    return re.sub(r"\s+", " ", normalize_text(text))


def clean_value(value: object) -> str:
    """Return a clean string; None / nan / none / null → empty string."""
    if value is None:
        return ""
    text = re.sub(r"\s+", " ", str(value)).strip()
    return "" if text.lower() in {"nan", "none", "null", ""} else text


def first(items: Iterable[str]) -> str:
    """Return the first non-empty clean string from an iterable."""
    for item in items:
        item = clean_value(item)
        if item:
            return item
    return ""


def dedupe(items: Iterable[str]) -> list[str]:
    """Return unique items preserving original order (case-insensitive key)."""
    seen: set[str] = set()
    output: list[str] = []
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
    """Convert a spelled-out number (e.g. 'twenty three') to an integer."""
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


# ===========================================================================
# DATA STRUCTURES
# ===========================================================================

@dataclass
class FieldResult:
    """Holds one extracted field's value plus quality metadata."""
    value: str
    confidence: float
    evidence: str   # verbatim snippet that led to this value
    source: str     # which part of the document was used


@dataclass
class DocumentContext:
    """All text views of a single PDF, built once and reused."""
    pdf_path: Path
    full_text: str
    first_page_text: str
    table_text: str
    sections: dict[str, str] = field(default_factory=dict)

    @property
    def combined_text(self) -> str:
        """Full text + table text joined."""
        return normalize_text(
            self.full_text + ("\n\n" + self.table_text if self.table_text else "")
        )

    @property
    def methods_text(self) -> str:
        """Methods + Results + Discussion concatenated."""
        return "\n".join(
            part for part in [
                self.sections.get("methods", ""),
                self.sections.get("results", ""),
                self.sections.get("discussion", ""),
            ] if part
        ).strip()


# ===========================================================================
# LOCAL LLM HELPER (optional Ollama back-end)
# ===========================================================================

class OllamaHelper:
    """
    Thin wrapper around the Ollama REST API for local LLM repair.
    Falls back silently if Ollama is not running or the model is unavailable.
    """

    def __init__(self, model: str, base_url: str = "http://127.0.0.1:11434"):
        self.model = model
        self.base_url = base_url.rstrip("/")

    def available(self) -> bool:
        """Return True only if Ollama is running and the model is loaded."""
        try:
            req = urllib.request.Request(f"{self.base_url}/api/tags", method="GET")
            with urllib.request.urlopen(req, timeout=3) as response:
                payload = json.loads(response.read().decode("utf-8"))
            return any(
                m.get("name", "").startswith(self.model)
                for m in payload.get("models", [])
            )
        except Exception:
            return False

    def ask(self, prompt: str, timeout: int = 120) -> str:
        """Send a prompt and return the model's response text."""
        body = json.dumps({
            "model": self.model,
            "prompt": prompt,
            "stream": False,
            "options": {"temperature": 0},
        }).encode("utf-8")
        req = urllib.request.Request(
            f"{self.base_url}/api/generate",
            data=body,
            headers={"Content-Type": "application/json"},
            method="POST",
        )
        with urllib.request.urlopen(req, timeout=timeout) as response:
            payload = json.loads(response.read().decode("utf-8"))
        return clean_value(payload.get("response", ""))


# ===========================================================================
# MAIN EXTRACTOR CLASS
# ===========================================================================

class GenomicsExtractor:
    """
    Extracts structured data from genome-wide gene family identification
    study PDFs using purely local, rule-based methods.

    Parameters
    ----------
    review_threshold : float
        Confidence below which a field is flagged for manual review.
        Range 0–1, default 0.65.
    use_tables : bool
        Enable Camelot-based table extraction (requires camelot-py[cv]).
    ollama_model : str or None
        Name of an Ollama model to use for repairing low-confidence fields.
        Pass None (default) to keep everything purely offline/rule-based.
    """

    def __init__(
        self,
        review_threshold: float = 0.65,
        use_tables: bool = False,
        ollama_model: Optional[str] = None,
    ) -> None:
        self.review_threshold = review_threshold
        self.use_tables = use_tables
        self.ollama = OllamaHelper(ollama_model) if ollama_model else None
        self.ollama_enabled = bool(self.ollama and self.ollama.available())

    # ------------------------------------------------------------------
    # PDF TEXT EXTRACTION
    # ------------------------------------------------------------------

    def _reader_score(self, text: str) -> float:
        """
        Score a text extraction by heuristics (word density, noise ratio,
        domain relevance). Used to choose the best reader when multiple
        back-ends are available.
        """
        lines = [line.strip() for line in text.splitlines() if line.strip()]
        if not lines:
            return float("-inf")
        sample = lines[:80]
        score: float = sum(len(re.findall(r"[A-Za-z]{2,}", ln)) for ln in sample)
        score -= sum(1 for ln in sample if len(ln) <= 2) * 5
        score -= sum(1 for ln in sample if re.fullmatch(r"[A-Za-z]", ln)) * 8
        # Reward domain-relevant signals
        if any("gene family" in ln.lower() for ln in lines[:25]):
            score += 15
        if any(ln.lower() == "abstract" for ln in lines[:35]):
            score += 10
        if any("genome-wide" in ln.lower() for ln in lines[:25]):
            score += 10
        return score

    def _extract_reader_outputs(self, pdf_path: Path) -> dict[str, str]:
        """Try each available PDF reader and return all successful outputs."""
        outputs: dict[str, str] = {}
        if pdfplumber is not None:
            try:
                with pdfplumber.open(str(pdf_path)) as pdf:
                    parts = [(page.extract_text() or "").strip() for page in pdf.pages]
                text = normalize_text("\n\n".join(p for p in parts if p))
                if text:
                    outputs["pdfplumber"] = text
            except Exception:
                pass
        if PdfReader is not None:
            try:
                reader = PdfReader(str(pdf_path))
                parts = [(page.extract_text() or "").strip() for page in reader.pages]
                text = normalize_text("\n\n".join(p for p in parts if p))
                if text:
                    outputs["pypdf"] = text
            except Exception:
                pass
        return outputs

    def _extract_best_text(self, pdf_path: Path) -> str:
        """Return the highest-scoring full-document text extraction."""
        outputs = self._extract_reader_outputs(pdf_path)
        return max(outputs.values(), key=self._reader_score) if outputs else ""

    def _extract_best_first_page(self, pdf_path: Path) -> str:
        """Return the highest-scoring first-page text extraction."""
        candidates: list[str] = []
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
        return max(candidates, key=self._reader_score) if candidates else ""

    def _extract_table_text(self, pdf_path: Path) -> str:
        """
        Extract table content using Camelot (optional).
        Returns empty string if Camelot is not installed or use_tables=False.
        """
        if not self.use_tables or camelot is None:
            return ""
        snippets: list[str] = []
        for flavor in ("lattice", "stream"):
            try:
                tables = camelot.read_pdf(str(pdf_path), pages="all", flavor=flavor)
            except Exception:
                continue
            for table in tables:
                try:
                    for row in table.df.fillna("").values.tolist():
                        joined = " | ".join(
                            str(cell).strip() for cell in row if str(cell).strip()
                        )
                        if joined:
                            snippets.append(joined)
                except Exception:
                    continue
            if snippets:
                break
        return normalize_text("\n".join(snippets))

    # ------------------------------------------------------------------
    # DOCUMENT STRUCTURE
    # ------------------------------------------------------------------

    def _lines(self, text: str, limit: Optional[int] = None) -> list[str]:
        """Return non-empty, stripped lines from text."""
        result = [ln.strip() for ln in normalize_text(text).splitlines() if ln.strip()]
        return result[:limit] if limit else result

    def _split_sections(self, text: str) -> dict[str, str]:
        """
        Split the document into named sections by recognising standard
        academic-paper headings.
        """
        found: dict[str, list[str]] = {name: [] for name in SECTION_HEADERS}
        current: Optional[str] = None
        for line in self._lines(text):
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
        return {name: "\n".join(lines).strip() for name, lines in found.items()}

    def build_context(self, pdf_path: Path) -> DocumentContext:
        """Extract all text views from a PDF and return a DocumentContext."""
        full_text = self._extract_best_text(pdf_path)
        first_page = self._extract_best_first_page(pdf_path) or full_text[:8000]
        table_text = self._extract_table_text(pdf_path)
        combined = normalize_text(full_text + ("\n\n" + table_text if table_text else ""))
        return DocumentContext(
            pdf_path=pdf_path,
            full_text=full_text,
            first_page_text=first_page,
            table_text=table_text,
            sections=self._split_sections(combined),
        )

    # ------------------------------------------------------------------
    # EVIDENCE EXTRACTION HELPERS
    # ------------------------------------------------------------------

    def _evidence(self, text: str, patterns: Iterable[str]) -> str:
        """Return the first matching line snippet (max 220 chars)."""
        for line in self._lines(text, 300):
            if any(pat.lower() in line.lower() for pat in patterns):
                return clean_value(line)[:220]
        return clean_value(first(self._lines(text, 20)))[:220]

    def _looks_like_journal(self, line: str) -> bool:
        low = line.lower()
        return (
            line in JOURNAL_HINTS
            or "journal homepage" in low
            or re.search(
                r"\b(volume|issn|eissn|research article|original article)\b", low
            ) is not None
            or any(hint.lower() in low for hint in JOURNAL_HINTS)
        )

    def _looks_like_author_line(self, line: str) -> bool:
        low = line.lower()
        if any(token in low for token in NOISE_TOKENS):
            return False
        if any(
            token in low
            for token in ["university", "department", "college", "institute", "faculty"]
        ):
            return False
        if "@" in line or self._looks_like_journal(line):
            return False
        name_hits = len(re.findall(r"[A-Z][a-z]+(?:\s+[A-Z][a-z]+)+", line))
        return name_hits >= 2 or (
            ";" in line and any(ch.isalpha() for ch in line)
        )

    @staticmethod
    def _clean_author_line(line: str) -> str:
        """Strip superscript numerals and noise characters from author lines."""
        line = re.sub(r"(?<=[A-Za-z])\d+", "", line)   # footnote numbers
        line = re.sub(r"[*#\u2020†‡§¶]+", "", line)    # symbols
        line = re.sub(r"\s+", " ", line)
        return line.strip(" ,;")

    # ------------------------------------------------------------------
    # FIELD-SPECIFIC EXTRACTORS
    # ------------------------------------------------------------------

    def extract_title(self, first_page_text: str) -> str:
        """Identify the paper title from the first page."""
        candidates: list[tuple[int, str, int]] = []
        front = self._lines(first_page_text, 40)
        for idx, line in enumerate(front):
            low = line.lower()
            if len(line) < 20 or len(line) > 260:
                continue
            if any(token in low for token in NOISE_TOKENS):
                continue
            if self._looks_like_author_line(line) or self._looks_like_journal(line):
                continue
            score = 0
            if idx <= 10:
                score += 5
            if re.search(
                r"\b(genome-wide|gene family|identification|analysis|"
                r"characterization|evolutionary|phylogenetic|expression)\b", low
            ):
                score += 8
            if re.search(r"\b(transcription factor|protein kinase|receptor)\b", low):
                score += 4
            if line[:1].isupper():
                score += 2
            if re.search(r"[.!?]$", line):
                score -= 2
            if re.search(r"\b(volume|issue|issn|eissn)\b", low):
                score -= 6
            candidates.append((score, line, idx))

        if not candidates:
            return ""
        candidates.sort(key=lambda item: item[0], reverse=True)
        best_line = candidates[0][1]
        best_idx = candidates[0][2]

        # Try to merge a continuation line (e.g. title wrapped to two lines)
        merged = [best_line]
        for line in front[best_idx + 1: best_idx + 5]:
            if (
                self._looks_like_author_line(line)
                or self._looks_like_journal(line)
                or any(token in line.lower() for token in NOISE_TOKENS)
                or len(line) > 130
            ):
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
        """Extract the author line from the first page."""
        front = self._lines(first_page_text, 45)
        for idx, line in enumerate(front):
            if title and clean_value(line) == clean_value(title):
                for follow in front[idx + 1: idx + 10]:
                    if self._looks_like_author_line(follow):
                        return self._clean_author_line(follow)
            if self._looks_like_author_line(line):
                return self._clean_author_line(line)
        return ""

    def extract_year(self, text: str) -> str:
        """Return the most frequent plausible publication year."""
        years = re.findall(r"\b(19\d{2}|20\d{2})\b", text[:12000])
        valid = [y for y in years if 1990 <= int(y) <= 2035]
        return Counter(valid).most_common(1)[0][0] if valid else ""

    def extract_journal(self, first_page_text: str) -> str:
        """Identify the journal name from the first page."""
        front = self._lines(first_page_text, 30)
        # Exact match against known journal list
        for line in front:
            if line in JOURNAL_HINTS:
                return line
        # Pattern: "In <JournalName> (year)" or "JournalName. 2023"
        for line in front:
            match = re.search(r"\.([A-Z][A-Za-z& \-]+)\s+(19|20)\d{2}\b", line)
            if match:
                return clean_value(match.group(1))
        # Structural hint
        for line in front:
            if (
                self._looks_like_journal(line)
                and line.lower() not in {"article", "research article", "original article"}
            ):
                return clean_value(line)
        return ""

    def extract_doi(self, text: str) -> str:
        """Extract DOI using the standard 10.XXXX/... pattern."""
        match = re.search(
            r"\b10\.\d{4,9}/[-._;()/:A-Z0-9]+\b", text, flags=re.IGNORECASE
        )
        return clean_value(match.group(0)) if match else ""

    def extract_keywords(self, ctx: DocumentContext) -> str:
        """Extract author-supplied keywords."""
        keywords = ctx.sections.get("keywords", "")
        if not keywords:
            match = re.search(
                r"keywords?\s*[:\-]\s*(.{5,300})", ctx.combined_text,
                flags=re.IGNORECASE
            )
            keywords = match.group(1) if match else ""
        # Keep only the first line to avoid absorbing body text
        keywords = re.split(r"\n|\. ", keywords)[0]
        if (
            "identified" in keywords.lower()
            or len(re.findall(r"[A-Za-z]{4,}", keywords)) > 20
        ):
            return ""
        return ", ".join(dedupe(re.split(r"[;,]", keywords)[:15]))

    def extract_country(self, first_page_text: str) -> str:
        """Identify the country of the corresponding/first author."""
        found: list[str] = []
        for country in COUNTRIES:
            if re.search(
                r"\b" + re.escape(country) + r"\b",
                first_page_text, flags=re.IGNORECASE
            ):
                found.append("United States" if country == "USA" else country)
        return Counter(found).most_common(1)[0][0] if found else ""

    def extract_organisms(
        self, title: str, abstract: str, keywords: str
    ) -> tuple[str, str]:
        """
        Identify primary and secondary study organisms.

        Strategy:
        1. Exact match against COMMON_SPECIES (weighted ×3).
        2. Regex for binomial names (Genus species).
        3. Return top-2 by frequency; filter noise.
        """
        context = "\n".join([title, abstract, keywords])
        found: list[str] = []

        # High-weight exact hits from known list
        for species in COMMON_SPECIES:
            if species.lower() in context.lower():
                found.extend([species] * 3)

        # General binomial name pattern
        found.extend(re.findall(r"\b([A-Z][a-z]+\s+[a-z][a-z\-]+)\b", context))

        _NOISE_STARTS = (
            "this study", "many studies", "overexpression of",
            "these results", "the present", "our results",
        )
        filtered = [
            item for item in found
            if len(item.split()) == 2
            and not item.lower().startswith(_NOISE_STARTS)
        ]
        common = [name for name, _ in Counter(filtered).most_common(2)]
        return (
            common[0] if common else "",
            common[1] if len(common) > 1 else "",
        )

    def extract_gene_family(self, title: str, abstract: str) -> str:
        """Extract the gene family name from title and abstract."""
        search_space = title + "\n" + abstract
        # Pattern-based extraction
        _PATTERNS = [
            r"\b([A-Za-z0-9/\-]+)\s+gene family\b",
            r"\b([A-Za-z0-9/\-]+)\s+transcription factors?\b",
            r"\b([A-Za-z0-9/\-]+)\s+family genes?\b",
            r"\b([A-Za-z0-9/\-]+)\s+(?:like\s+)?gene family\b",
            r"\b([A-Za-z0-9/\-]+)\s+genes?\s+in\s+[A-Z][a-z]+\b",
            r"genome-wide\s+(?:analysis|identification|characterization)\s+of\s+"
            r"([A-Za-z0-9/\-]+)\b",
        ]
        for pattern in _PATTERNS:
            match = re.search(pattern, search_space, flags=re.IGNORECASE)
            if match:
                value = re.sub(
                    r"^(the|a|an)\s+", "", match.group(1), flags=re.IGNORECASE
                ).strip()
                if value and len(value) >= 2:
                    return clean_value(value)

        # Fallback: known family name in title/abstract
        for family in KNOWN_GENE_FAMILIES:
            if re.search(
                r"\b" + re.escape(family) + r"\b",
                search_space, flags=re.IGNORECASE
            ):
                return family
        return ""

    def extract_genome_source(self, text: str) -> str:
        """Identify genome databases used as sequence sources."""
        _SOURCES = [
            "Phytozome", "NCBI", "Ensembl", "EnsemblPlants", "TAIR",
            "Gramene", "JGI", "MaizeGDB", "SGN", "RAP-DB", "CoGe",
            "Sol Genomics", "PLAZA", "CottonFGD",
        ]
        return ", ".join(
            dedupe(src for src in _SOURCES if src.lower() in text.lower())
        )

    def extract_hmmer(self, text: str) -> str:
        """Return HMMER usage + Pfam domain accessions if found."""
        used = re.search(
            r"\b(HMMER|hmmsearch|hmmscan|hmmprofile)\b",
            text, flags=re.IGNORECASE
        ) is not None
        pfams = dedupe(re.findall(r"\bPF\d{5}\b", text))
        has_pfam = re.search(r"\bPfam\b", text, flags=re.IGNORECASE) is not None
        if not used and not pfams and not has_pfam:
            return "No"
        if pfams:
            return "Yes; Pfam domains: " + ", ".join(pfams[:20])
        if has_pfam:
            return "Yes; Pfam mentioned"
        return "Yes"

    def extract_blast(self, text: str) -> str:
        """Return BLAST E-value cut-off or 'Yes'/'No'."""
        if re.search(r"\bBLAST[PNX]?\b", text, flags=re.IGNORECASE) is None:
            return "No"
        flat = flatten(text)
        _EVALUE_PATTERNS = [
            r"BLAST[PNX]?\b.{0,120}?E-?value(?:\s+cut-?off|\s+threshold)?"
            r"(?:\s+of)?\s*[:=<>~]?\s*([0-9.]+(?:e-?\d+)?)",
            r"E-?value(?:\s+cut-?off|\s+threshold)?(?:\s+of)?\s*[:=<>~]?\s*"
            r"([0-9.]+(?:e-?\d+)?)",
            r"([0-9.]+\s*[×x]\s*10[-\u2212]?\d+)",
        ]
        for pattern in _EVALUE_PATTERNS:
            match = re.search(pattern, flat, flags=re.IGNORECASE)
            if match:
                val = (
                    match.group(1)
                    .replace(" ", "")
                    .replace("×", "e")
                    .replace("x", "e")
                    .replace("10-", "e-")
                )
                return val
        return "Yes"

    def _tool_yes_no(self, text: str, aliases: Iterable[str]) -> str:
        """Return 'Yes'/'No' based on whether any alias appears in text."""
        return (
            "Yes"
            if any(
                re.search(r"\b" + re.escape(alias) + r"\b", text, flags=re.IGNORECASE)
                for alias in aliases
            )
            else "No"
        )

    def _extract_other_tools(
        self, text: str, pool: Iterable[str], exclude: Iterable[str]
    ) -> str:
        """Return comma-separated pool tools found in text, minus excluded ones."""
        exclude_set = {item.lower() for item in exclude}
        return ", ".join(
            dedupe(
                tool for tool in pool
                if tool.lower() not in exclude_set
                and re.search(
                    r"\b" + re.escape(tool) + r"\b", text, flags=re.IGNORECASE
                )
            )
        )

    def extract_qpcr(self, text: str) -> str:
        """Detect experimental qPCR validation."""
        return (
            "Yes"
            if re.search(
                r"\b(qPCR|RT-qPCR|qRT-PCR|real-time PCR|SYBR\s*Green|"
                r"quantitative\s*RT-PCR)\b",
                text, flags=re.IGNORECASE,
            )
            else "No"
        )

    def extract_gene_count(
        self, title: str, abstract: str, results: str
    ) -> str:
        """Extract the total number of gene family members identified."""
        search_space = flatten("\n".join([title, abstract, results]))
        _PATTERNS = [
            r"\ba total of ((?:\d{1,4})|(?:[a-z\-]+))\s+\S*\s*genes?\s+(?:were\s+)?identified\b",
            r"\bidentified\s+((?:\d{1,4})|(?:[a-z\-]+))\s+\S*\s*genes?\b",
            r"\bfound\s+((?:\d{1,4})|(?:[a-z\-]+))\s+\S*\s*genes?\b",
            r"\b((?:\d{1,4})|(?:[a-z\-]+))\s+members?\s+of\s+the\s+\S+\s+gene family\b",
            r"\bcomprises?\s+((?:\d{1,4})|(?:[a-z\-]+))\s+\S*\s*genes?\b",
            r"\bcontains?\s+((?:\d{1,4})|(?:[a-z\-]+))\s+\S*\s*genes?\b",
            r"\bconsists?\s+of\s+((?:\d{1,4})|(?:[a-z\-]+))\s+\S*\s*genes?\b",
        ]
        for pattern in _PATTERNS:
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
        """Extract the number and type of phylogenetic groupings."""
        flat = flatten(text)
        for pattern in [
            r"\bclassified\s+into\s+(\d+)\s+(clades?|groups?|subfamilies|clusters?)\b",
            r"\bdivided\s+into\s+(\d+)\s+(clades?|groups?|subfamilies|clusters?)\b",
            r"\b(\d+)\s+(clades?|groups?|subfamilies|clusters?)\b",
        ]:
            match = re.search(pattern, flat, flags=re.IGNORECASE)
            if match:
                return f"{match.group(1)} {match.group(2)}"
        return ""

    def extract_duplication(self, text: str) -> str:
        """Detect tandem and/or segmental duplication events."""
        events: list[str] = []
        if re.search(r"\btandem duplication\b", text, flags=re.IGNORECASE):
            events.append("TANDEM")
        if re.search(r"\bsegmental duplication\b", text, flags=re.IGNORECASE):
            events.append("SEGMENTAL")
        if re.search(r"\bwhole-genome duplication\b|WGD", text, flags=re.IGNORECASE):
            events.append("WGD")
        return "/".join(events)

    def extract_rnaseq(self, text: str) -> str:
        """Detect RNA-seq / transcriptome expression analysis."""
        return (
            "Yes"
            if re.search(
                r"\b(RNA-seq|RNA\s*sequencing|transcriptome|FPKM|TPM|RPKM|"
                r"DESeq|edgeR|HTSeq|featureCounts)\b",
                text, flags=re.IGNORECASE,
            )
            else "No"
        )

    # ------------------------------------------------------------------
    # CONFIDENCE SCORING
    # ------------------------------------------------------------------

    def _confidence(self, field: str, value: str) -> float:
        """
        Heuristic confidence score in [0, 1] for an extracted value.
        Higher = more reliable; fields below review_threshold are flagged.
        """
        value = clean_value(value)
        if not value:
            return 0.0

        if field == "DOI":
            return 0.98 if re.match(r"^10\.\d{4,9}/", value) else 0.30

        if field == "Year":
            return 0.95 if re.fullmatch(r"(19|20)\d{2}", value) else 0.25

        if field in {
            "SMART", "CDD", "MEGA", "IQ-TREE", "MEME", "MCScanX", "TBTools",
            "Experimental Validation (qPCR)", "Expression Analysis (RNA-seq)",
        }:
            return 0.90 if value in {"Yes", "No"} else 0.40

        if field == "Number of Genes Identified":
            return 0.92 if re.fullmatch(r"\d{1,4}", value) else 0.30

        if field == "HMMER (with Pfam domains)":
            if "PF" in value:
                return 0.95
            return 0.85 if "Yes" in value else 0.90

        if field == "BLAST (E-value cut-off used)":
            return 0.92 if re.search(r"e-?\d+", value, re.IGNORECASE) else (
                0.85 if value == "Yes" else 0.90
            )

        score = 0.55
        if len(value) >= 5:
            score += 0.10
        if field in {"Title", "Authors", "Journal"}:
            score += 0.10
        if re.search(
            r"\b(revised|published|citation|license|editor)\b",
            value, flags=re.IGNORECASE
        ):
            score -= 0.45
        if field in {"Primary Organism", "Secondary Organism"}:
            if len(value.split()) != 2:
                score -= 0.20
            elif value in COMMON_SPECIES:
                score += 0.15
        if field == "Journal" and re.search(r"\bissn\b", value, flags=re.IGNORECASE):
            score -= 0.30
        if field == "Gene Family" and value in KNOWN_GENE_FAMILIES:
            score += 0.15
        return round(max(0.0, min(1.0, score)), 2)

    # ------------------------------------------------------------------
    # RESULT FACTORY
    # ------------------------------------------------------------------

    def _make_result(
        self,
        field: str,
        value: str,
        evidence_text: str,
        source: str,
        patterns: Iterable[str],
    ) -> FieldResult:
        return FieldResult(
            value=clean_value(value),
            confidence=self._confidence(field, value),
            evidence=self._evidence(evidence_text, patterns),
            source=source,
        )

    # ------------------------------------------------------------------
    # OPTIONAL LOCAL LLM REPAIR
    # ------------------------------------------------------------------

    def _llm_repair(
        self, field: str, ctx: DocumentContext, current: FieldResult
    ) -> FieldResult:
        """
        Use Ollama to repair a low-confidence field.
        Returns the current result unchanged if Ollama is unavailable
        or if the current confidence already meets the threshold.
        """
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
            f"You are a bioinformatics data extractor.\n"
            f"Extract the value for the field '{field}' from the paper text below.\n"
            f"Rules:\n"
            f"- Return ONLY the value; no explanation.\n"
            f"- If the field is not present, return exactly: EMPTY\n"
            f"- For Yes/No fields, return only 'Yes' or 'No'.\n\n"
            f"Paper text:\n{evidence_block}"
        )
        try:
            answer = self.ollama.ask(prompt)  # type: ignore[union-attr]
        except (urllib.error.URLError, TimeoutError, OSError, ValueError):
            return current
        if not answer or answer.upper() == "EMPTY":
            return current
        repaired = self._make_result(field, answer, evidence_block, "ollama", [field])
        return repaired if repaired.confidence > current.confidence else current

    # ------------------------------------------------------------------
    # MAIN EXTRACTION PIPELINE
    # ------------------------------------------------------------------

    def extract_fields(self, ctx: DocumentContext) -> dict[str, FieldResult]:
        """
        Run all field extractors on the given DocumentContext.
        Returns a dict mapping column name → FieldResult.
        """
        title      = self.extract_title(ctx.first_page_text)
        authors    = self.extract_authors(ctx.first_page_text, title)
        year       = self.extract_year(ctx.first_page_text or ctx.full_text)
        journal    = self.extract_journal(ctx.first_page_text)
        doi        = self.extract_doi(ctx.combined_text)
        keywords   = self.extract_keywords(ctx)
        abstract   = ctx.sections.get("abstract", "")
        primary, secondary = self.extract_organisms(title, abstract, keywords)
        gene_family   = self.extract_gene_family(title, abstract)
        genome_source = self.extract_genome_source(ctx.combined_text)
        methods_text  = ctx.methods_text or ctx.combined_text

        data: dict[str, FieldResult] = {
            "Country of Origin": self._make_result(
                "Country of Origin",
                self.extract_country(ctx.first_page_text),
                ctx.first_page_text, "front", COUNTRIES,
            ),
            "Title": self._make_result(
                "Title", title, ctx.first_page_text, "front", [title]
            ),
            "Authors": self._make_result(
                "Authors", authors, ctx.first_page_text, "front", [authors]
            ),
            "Year": self._make_result(
                "Year", year, ctx.first_page_text, "front",
                ["received", "published", "accepted"],
            ),
            "Journal": self._make_result(
                "Journal", journal, ctx.first_page_text, "front", JOURNAL_HINTS
            ),
            "DOI": self._make_result(
                "DOI", doi, ctx.combined_text, "full",
                ["doi", "https://doi.org"],
            ),
            "Keywords": self._make_result(
                "Keywords", keywords,
                ctx.sections.get("keywords", "") or ctx.combined_text,
                "keywords", ["keywords"],
            ),
            "Primary Organism": self._make_result(
                "Primary Organism", primary,
                "\n".join([title, abstract, keywords]),
                "context", COMMON_SPECIES,
            ),
            "Secondary Organism": self._make_result(
                "Secondary Organism", secondary,
                "\n".join([title, abstract, keywords]),
                "context", COMMON_SPECIES,
            ),
            "Gene Family": self._make_result(
                "Gene Family", gene_family,
                "\n".join([title, abstract]),
                "front+abstract", KNOWN_GENE_FAMILIES,
            ),
            "Genome Source": self._make_result(
                "Genome Source", genome_source, ctx.combined_text, "full",
                ["Phytozome", "NCBI", "Ensembl", "TAIR", "JGI"],
            ),
            "HMMER (with Pfam domains)": self._make_result(
                "HMMER (with Pfam domains)",
                self.extract_hmmer(ctx.combined_text),
                ctx.combined_text, "full", ["HMMER", "Pfam"],
            ),
            "BLAST (E-value cut-off used)": self._make_result(
                "BLAST (E-value cut-off used)",
                self.extract_blast(ctx.combined_text),
                ctx.combined_text, "full", ["BLAST", "E-value"],
            ),
            "SMART": self._make_result(
                "SMART", self._tool_yes_no(ctx.combined_text, ["SMART"]),
                ctx.combined_text, "full", ["SMART"],
            ),
            "Other (Tools for Domain Analysis)": self._make_result(
                "Other (Tools for Domain Analysis)",
                self._extract_other_tools(
                    ctx.combined_text, DOMAIN_TOOLS, ["HMMER", "Pfam", "SMART", "CDD"]
                ),
                ctx.combined_text, "full", DOMAIN_TOOLS,
            ),
            "CDD": self._make_result(
                "CDD",
                self._tool_yes_no(
                    ctx.combined_text,
                    ["CDD", "Conserved Domain Database", "NCBI-CDD"],
                ),
                ctx.combined_text, "full", ["CDD"],
            ),
            "MEGA": self._make_result(
                "MEGA", self._tool_yes_no(ctx.combined_text, ["MEGA"]),
                ctx.combined_text, "full", ["MEGA"],
            ),
            "IQ-TREE": self._make_result(
                "IQ-TREE",
                self._tool_yes_no(
                    ctx.combined_text, ["IQ-TREE", "IQTREE", "iqtree"]
                ),
                ctx.combined_text, "full", ["IQ-TREE"],
            ),
            "MEME": self._make_result(
                "MEME", self._tool_yes_no(ctx.combined_text, ["MEME"]),
                ctx.combined_text, "full", ["MEME"],
            ),
            "Others (Tools for Motif Analysis)": self._make_result(
                "Others (Tools for Motif Analysis)",
                self._extract_other_tools(ctx.combined_text, MOTIF_TOOLS, ["MEME"]),
                ctx.combined_text, "full", MOTIF_TOOLS,
            ),
            "MCScanX": self._make_result(
                "MCScanX", self._tool_yes_no(ctx.combined_text, ["MCScanX"]),
                ctx.combined_text, "full", ["MCScanX"],
            ),
            "Other (Tools for Synteny)": self._make_result(
                "Other (Tools for Synteny)",
                self._extract_other_tools(
                    ctx.combined_text, SYNTENY_TOOLS, ["MCScanX"]
                ),
                ctx.combined_text, "full", SYNTENY_TOOLS,
            ),
            "TBTools": self._make_result(
                "TBTools",
                self._tool_yes_no(ctx.combined_text, ["TBtools", "TBTools"]),
                ctx.combined_text, "full", ["TBtools"],
            ),
            "Other (Tools for Visualization)": self._make_result(
                "Other (Tools for Visualization)",
                self._extract_other_tools(ctx.combined_text, VIS_TOOLS, ["TBtools"]),
                ctx.combined_text, "full", VIS_TOOLS,
            ),
            "Experimental Validation (qPCR)": self._make_result(
                "Experimental Validation (qPCR)",
                self.extract_qpcr(methods_text),
                methods_text, "methods", ["qPCR", "RT-qPCR"],
            ),
            "Other Tools Used In Methodology": self._make_result(
                "Other Tools Used In Methodology",
                self._extract_other_tools(methods_text, METHOD_TOOLS, []),
                methods_text, "methods", METHOD_TOOLS,
            ),
            "Number of Genes Identified": self._make_result(
                "Number of Genes Identified",
                self.extract_gene_count(
                    title, abstract, ctx.sections.get("results", "")
                ),
                "\n".join([title, abstract, ctx.sections.get("results", "")]),
                "title+abstract+results",
                ["identified", "genes"],
            ),
            "Phylogenetic Structure": self._make_result(
                "Phylogenetic Structure",
                self.extract_phylo_structure(ctx.combined_text),
                ctx.combined_text, "full",
                ["clades", "groups", "subfamilies"],
            ),
            "Gene Duplication Events (TANDEM/SEGMENTAL)": self._make_result(
                "Gene Duplication Events (TANDEM/SEGMENTAL)",
                self.extract_duplication(ctx.combined_text),
                ctx.combined_text, "full",
                ["tandem duplication", "segmental duplication"],
            ),
            "Expression Analysis (RNA-seq)": self._make_result(
                "Expression Analysis (RNA-seq)",
                self.extract_rnaseq(methods_text),
                methods_text, "methods",
                ["RNA-seq", "transcriptome"],
            ),
        }

        # Optional LLM repair pass for key metadata fields
        if self.ollama_enabled:
            _LLM_REPAIR_FIELDS = [
                "Title", "Authors", "Journal", "Primary Organism",
                "Secondary Organism", "Gene Family", "Number of Genes Identified",
                "Country of Origin", "Year",
            ]
            for f in _LLM_REPAIR_FIELDS:
                data[f] = self._llm_repair(f, ctx, data[f])

        return data

    def _flatten_row(self, fields: dict[str, FieldResult]) -> dict[str, str]:
        """Return a plain dict of column → value string."""
        return {col: fields[col].value for col in COLUMNS}

    def _review_rows(
        self, pdf_name: str, fields: dict[str, FieldResult]
    ) -> list[dict]:
        """Return rows for the review queue CSV (low-confidence fields only)."""
        return [
            {
                "paper": pdf_name,
                "field": field,
                "value": result.value,
                "confidence": result.confidence,
                "source": result.source,
                "evidence": result.evidence,
            }
            for field, result in fields.items()
            if result.confidence < self.review_threshold
        ]

    def _summary_row(
        self, pdf_name: str, fields: dict[str, FieldResult]
    ) -> dict:
        """Return the per-paper quality summary row."""
        confidences = [r.confidence for r in fields.values()]
        low = [f for f, r in fields.items() if r.confidence < self.review_threshold]
        return {
            "paper": pdf_name,
            "average_confidence": round(sum(confidences) / len(confidences), 2),
            "low_confidence_fields": ", ".join(low),
            "needs_review": "Yes" if low else "No",
        }

    # ------------------------------------------------------------------
    # BATCH PROCESSING
    # ------------------------------------------------------------------

    def process_directory(
        self,
        input_dir: Path,
        output_dir: Path,
        run_name: str,
        excel: bool,
    ) -> None:
        """
        Process all PDFs in input_dir and write results to output_dir.

        Output files
        ------------
        <run_name>.csv              Main extraction table.
        <run_name>_review_queue.csv Low-confidence cells for manual review.
        <run_name>_paper_summary.csv Per-paper quality summary.
        <run_name>_evidence.jsonl   Field-level evidence for traceability.
        <run_name>_README.txt       Run guide.
        <run_name>.xlsx             (optional) Excel copy of the main CSV.
        """
        pdf_files = sorted(input_dir.rglob("*.pdf"))
        if not pdf_files:
            raise FileNotFoundError(f"No PDF files found in {input_dir}")

        output_dir.mkdir(parents=True, exist_ok=True)

        main_csv      = output_dir / f"{run_name}.csv"
        review_csv    = output_dir / f"{run_name}_review_queue.csv"
        summary_csv   = output_dir / f"{run_name}_paper_summary.csv"
        evidence_jsonl = output_dir / f"{run_name}_evidence.jsonl"
        readme        = output_dir / f"{run_name}_README.txt"

        rows: list[dict]   = []
        reviews: list[dict] = []
        summaries: list[dict] = []
        evidence_lines: list[str] = []

        print(f"\nGenomics PDF Extractor")
        print(f"======================")
        print(f"Input folder : {input_dir}")
        print(f"Output folder: {output_dir}")
        print(f"PDFs found   : {len(pdf_files)}")
        ollama_label = f"Yes ({self.ollama.model})" if self.ollama_enabled else "No"
        print(f"Ollama repair: {ollama_label}")
        print(f"Table extract: {'Yes' if self.use_tables else 'No'}\n")

        for idx, pdf_path in enumerate(pdf_files, start=1):
            print(f"[{idx:>3}/{len(pdf_files)}] {pdf_path.name}")
            try:
                ctx    = self.build_context(pdf_path)
                fields = self.extract_fields(ctx)
                rows.append(self._flatten_row(fields))
                reviews.extend(self._review_rows(pdf_path.name, fields))
                summaries.append(self._summary_row(pdf_path.name, fields))
                evidence_lines.append(
                    json.dumps(
                        {
                            "paper": pdf_path.name,
                            "fields": {
                                f: {
                                    "value":      r.value,
                                    "confidence": r.confidence,
                                    "evidence":   r.evidence,
                                    "source":     r.source,
                                }
                                for f, r in fields.items()
                            },
                        },
                        ensure_ascii=False,
                    )
                )
                # Print gene family + count as quick QC
                gf  = fields["Gene Family"].value
                cnt = fields["Number of Genes Identified"].value
                org = fields["Primary Organism"].value
                print(f"         → {org or '?'} | {gf or '?'} | {cnt or '?'} genes")
            except Exception as exc:
                print(f"         ✗ FAILED: {exc}")
                rows.append({col: "" for col in COLUMNS})
                summaries.append({
                    "paper": pdf_path.name,
                    "average_confidence": 0.0,
                    "low_confidence_fields": "ALL_FIELDS",
                    "needs_review": "Yes",
                })

        # Write outputs
        df = pd.DataFrame(rows, columns=COLUMNS)
        df.to_csv(main_csv, index=False, encoding="utf-8-sig")

        pd.DataFrame(
            reviews,
            columns=["paper", "field", "value", "confidence", "source", "evidence"],
        ).to_csv(review_csv, index=False, encoding="utf-8-sig")

        pd.DataFrame(
            summaries,
            columns=[
                "paper", "average_confidence",
                "low_confidence_fields", "needs_review",
            ],
        ).to_csv(summary_csv, index=False, encoding="utf-8-sig")

        evidence_jsonl.write_text("\n".join(evidence_lines), encoding="utf-8")

        readme.write_text(
            "\n".join([
                "Genomics PDF Extractor — Run Report",
                "====================================",
                f"Input folder   : {input_dir}",
                f"Total PDFs     : {len(pdf_files)}",
                f"Main CSV       : {main_csv.name}",
                f"Review queue   : {review_csv.name}",
                f"Paper summary  : {summary_csv.name}",
                f"Evidence JSONL : {evidence_jsonl.name}",
                f"Review threshold: {self.review_threshold}",
                f"Table extraction: {'Yes' if self.use_tables else 'No'}",
                f"Ollama model   : {self.ollama.model if self.ollama_enabled else 'None'}",
                "",
                "Recommended workflow:",
                "1. Open the main CSV for first-pass data.",
                "2. Check paper_summary.csv for average confidence per paper.",
                "3. Open review_queue.csv to manually verify low-confidence cells.",
                "4. Use evidence.jsonl for traceable field-by-field provenance.",
            ]),
            encoding="utf-8",
        )

        if excel:
            xlsx_path = output_dir / f"{run_name}.xlsx"
            df.to_excel(xlsx_path, index=False)
            print(f"\nSaved Excel       : {xlsx_path}")

        print(f"\nSaved CSV         : {main_csv}")
        print(f"Saved review queue: {review_csv}")
        print(f"Saved paper summary: {summary_csv}")
        print(f"Saved evidence    : {evidence_jsonl}")
        print(f"Saved run guide   : {readme}")


    def process_single_file(
        self,
        pdf_path: Path,
        output_dir: Path,
        run_name: str,
        excel: bool,
        verbose: bool = False,
    ) -> None:
        """
        Process a single PDF and write results to output_dir.
        Identical output structure to process_directory but for one file.
        """
        if not pdf_path.exists():
            raise FileNotFoundError(f"PDF not found: {pdf_path}")
        if pdf_path.suffix.lower() != ".pdf":
            raise ValueError(f"Not a PDF file: {pdf_path}")

        output_dir.mkdir(parents=True, exist_ok=True)

        print("\nGenomics PDF Extractor (single-file mode)")
        print("==========================================")
        print(f"File         : {pdf_path}")
        print(f"Output folder: {output_dir}")
        ollama_label = f"Yes ({self.ollama.model})" if self.ollama_enabled else "No"
        print(f"Ollama repair: {ollama_label}")
        table_label = "Yes" if self.use_tables else "No"
        print(f"Table extract: {table_label}\n")

        ctx    = self.build_context(pdf_path)
        fields = self.extract_fields(ctx)

        if verbose:
            print("\n── Extraction results ────────────────────────────────")
            for col in COLUMNS:
                r = fields[col]
                flag = " ⚠" if r.confidence < self.review_threshold else ""
                print(f"  {col:<45} conf={r.confidence:.2f}{flag}")
                print(f"    {r.value[:90]!r}")
            print()

        main_csv       = output_dir / f"{run_name}.csv"
        review_csv     = output_dir / f"{run_name}_review_queue.csv"
        summary_csv    = output_dir / f"{run_name}_paper_summary.csv"
        evidence_jsonl = output_dir / f"{run_name}_evidence.jsonl"

        import pandas as _pd
        _pd.DataFrame([self._flatten_row(fields)], columns=COLUMNS).to_csv(
            main_csv, index=False, encoding="utf-8-sig"
        )
        reviews = self._review_rows(pdf_path.name, fields)
        _pd.DataFrame(
            reviews,
            columns=["paper", "field", "value", "confidence", "source", "evidence"],
        ).to_csv(review_csv, index=False, encoding="utf-8-sig")
        _pd.DataFrame(
            [self._summary_row(pdf_path.name, fields)],
            columns=["paper", "average_confidence", "low_confidence_fields", "needs_review"],
        ).to_csv(summary_csv, index=False, encoding="utf-8-sig")
        evidence_jsonl.write_text(
            json.dumps({
                "paper": pdf_path.name,
                "fields": {
                    f: {"value": r.value, "confidence": r.confidence,
                        "evidence": r.evidence, "source": r.source}
                    for f, r in fields.items()
                },
            }, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
        if excel:
            xlsx_path = output_dir / f"{run_name}.xlsx"
            _pd.DataFrame([self._flatten_row(fields)], columns=COLUMNS).to_excel(
                xlsx_path, index=False
            )
            print(f"Saved Excel       : {xlsx_path}")
        print(f"Saved CSV         : {main_csv}")
        print(f"Saved review queue: {review_csv}")
        print(f"Saved paper summary: {summary_csv}")
        print(f"Saved evidence    : {evidence_jsonl}")

# ===========================================================================
# CLI
# ===========================================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Genomics PDF Extractor — fully local, no-API tool for extracting "
            "structured data from genome-wide gene family identification PDFs."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input", default=".",
        help="Folder containing PDF files (searched recursively).",
    )
    parser.add_argument(
        "--output-dir", default="extractor_output",
        help="Folder where all result files will be written.",
    )
    parser.add_argument(
        "--run-name", default="genomics_extraction",
        help="Base name used for all output files.",
    )
    parser.add_argument(
        "--excel", action="store_true",
        help="Also save an Excel (.xlsx) copy of the main CSV.",
    )
    parser.add_argument(
        "--review-threshold", type=float, default=0.65,
        help="Confidence cut-off below which a field is flagged for review.",
    )
    parser.add_argument(
        "--use-tables", action="store_true",
        help="Enable Camelot table extraction (requires camelot-py[cv]).",
    )
    parser.add_argument(
        "--ollama-model", default="",
        help=(
            "Optional local Ollama model for low-confidence field repair. "
            "Example: qwen2.5:7b  |  Requires Ollama running locally."
        ),
    )
    parser.add_argument(
        "--file", default="",
        help=(
            "Process a single PDF instead of a whole directory. "
            "Mutually exclusive with --input."
        ),
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="Print every extracted field and confidence score to stdout.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    ensure_dependencies()

    single_file = clean_value(args.file)
    if single_file and clean_value(args.input) not in {"", "."}:
        print("Error: --file and --input cannot be used together.", file=sys.stderr)
        return 1

    extractor = GenomicsExtractor(
        review_threshold=args.review_threshold,
        use_tables=args.use_tables,
        ollama_model=clean_value(args.ollama_model) or None,
    )
    run_name = clean_value(args.run_name) or "genomics_extraction"
    output_dir = Path(args.output_dir).expanduser().resolve()

    if single_file:
        extractor.process_single_file(
            pdf_path=Path(single_file).expanduser().resolve(),
            output_dir=output_dir,
            run_name=run_name,
            excel=args.excel,
            verbose=args.verbose,
        )
    else:
        extractor.process_directory(
            input_dir=Path(args.input).expanduser().resolve(),
            output_dir=output_dir,
            run_name=run_name,
            excel=args.excel,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
