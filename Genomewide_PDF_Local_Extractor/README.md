# genomewide_pdf_extractor

[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![No API Required](https://img.shields.io/badge/API-None%20required-green.svg)]()

A **fully local, zero-API** tool for extracting structured metadata from genome-wide gene family identification study PDFs. Designed for bioinformatics researchers performing systematic reviews or meta-analyses.

---

## Features

- **100% local** — no internet connection or API key required during extraction
- **Dual PDF reader** — tries both `pdfplumber` and `pypdf`, picks the best result automatically
- **30 structured fields** extracted per paper (title, authors, organisms, tools, gene counts, etc.)
- **Confidence scoring** — every field gets a quality score; low-confidence cells are flagged for review
- **Evidence tracing** — every extracted value is linked to the verbatim snippet it came from
- **Optional local LLM repair** — plug in any [Ollama](https://ollama.com) model to fix low-confidence fields without any cloud dependency
- **Optional table extraction** — Camelot back-end for structured table parsing
- **Batch processing** — drop a folder of PDFs and get a single structured CSV

---

## Extracted Fields

| Category | Fields |
|---|---|
| **Bibliographic** | Country of Origin, Title, Authors, Year, Journal, DOI, Keywords |
| **Biological** | Primary Organism, Secondary Organism, Gene Family, Genome Source |
| **Domain Analysis** | HMMER (with Pfam domains), BLAST (E-value), SMART, CDD, Other Domain Tools |
| **Phylogenetics** | MEGA, IQ-TREE, Phylogenetic Structure, Gene Duplication Events |
| **Motif Analysis** | MEME, Other Motif Tools |
| **Synteny** | MCScanX, Other Synteny Tools |
| **Visualization** | TBTools, Other Visualization Tools |
| **Experimental** | Experimental Validation (qPCR), Expression Analysis (RNA-seq) |
| **Results** | Number of Genes Identified, Other Tools Used |

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/YOUR_USERNAME/genomewide-pdf-extractor.git
cd genomewide-pdf-extractor
```

### 2. Install dependencies

```bash
pip install -r requirements.txt
```

### Optional: Camelot (table extraction)

Only needed if you plan to use `--use-tables`:

```bash
pip install camelot-py[cv]
# Also requires Ghostscript: https://www.ghostscript.com/download.html
```

### Optional: Ollama (local LLM repair)

Install [Ollama](https://ollama.com), then pull a model:

```bash
ollama pull qwen2.5:7b
```

---

## Quick Start

```bash
# Basic: extract from all PDFs in a folder
python genomewide_pdf_extractor.py --input ./my_papers --output-dir ./results

# Single file (useful for testing or spot-checking)
python genomewide_pdf_extractor.py --file ./my_papers/paper1.pdf --output-dir ./results

# Show all extracted fields and confidence scores in the terminal
python genomewide_pdf_extractor.py --file paper.pdf --verbose

# Also save an Excel file
python genomewide_pdf_extractor.py --input ./my_papers --output-dir ./results --excel

# With local LLM repair for low-confidence fields (Ollama required)
python genomewide_pdf_extractor.py --input ./my_papers --output-dir ./results \
    --ollama-model qwen2.5:7b

# Lower the review threshold (flag more fields for manual check)
python genomewide_pdf_extractor.py --input ./my_papers --review-threshold 0.75

# Enable table extraction (slower, sometimes more complete)
python genomewide_pdf_extractor.py --input ./my_papers --use-tables
```

---

## Output Files

| File | Description |
|---|---|
| `<run_name>.csv` | Main extraction table — one row per PDF |
| `<run_name>_review_queue.csv` | Low-confidence cells flagged for manual review |
| `<run_name>_paper_summary.csv` | Per-paper average confidence + which fields need review |
| `<run_name>_evidence.jsonl` | Full evidence trace — value + snippet + source per field |
| `<run_name>_README.txt` | Run guide generated at extraction time |
| `<run_name>.xlsx` | *(optional)* Excel copy of the main CSV |

### Understanding the review queue

The `_review_queue.csv` lists every field whose confidence score fell below `--review-threshold` (default: **0.65**). Open this file alongside the original PDF to quickly correct only the uncertain cells — you do not need to manually re-check high-confidence extractions.

---

## CLI Reference

```
usage: genomewide_pdf_extractor.py [-h] [--input INPUT] [--file FILE]
                              [--output-dir OUTPUT_DIR] [--run-name RUN_NAME]
                              [--excel] [--review-threshold THRESHOLD]
                              [--use-tables] [--ollama-model MODEL]
                              [--verbose]

options:
  --input               Folder with PDFs (searched recursively). Default: .
  --file                Process a single PDF file instead of a directory.
  --output-dir          Output folder. Default: extractor_output
  --run-name            Base name for output files. Default: genomics_extraction
  --excel               Also write an XLSX copy of the main CSV.
  --review-threshold    Confidence cut-off for review flagging. Default: 0.65
  --use-tables          Enable Camelot table extraction (requires camelot-py[cv]).
  --ollama-model        Ollama model name for LLM repair (e.g. qwen2.5:7b).
  --verbose             Print every field + confidence score to stdout.
```

---

## How It Works

```
PDF files
   │
   ▼
┌─────────────────────────────────┐
│  PDF Text Extraction            │
│  pdfplumber + pypdf (best wins) │
└────────────┬────────────────────┘
             │
             ▼
┌─────────────────────────────────┐
│  Document Structure Detection  │
│  Section splitter (Abstract,   │
│  Methods, Results, etc.)       │
└────────────┬────────────────────┘
             │
             ▼
┌─────────────────────────────────┐
│  Rule-Based Field Extraction   │
│  30 fields via regex + domain  │
│  vocabulary matching           │
└────────────┬────────────────────┘
             │
             ▼
┌─────────────────────────────────┐
│  Confidence Scoring            │
│  Per-field heuristic scoring   │
└────────────┬────────────────────┘
             │
   (optional)▼
┌─────────────────────────────────┐
│  Local LLM Repair (Ollama)     │
│  Only for low-confidence fields│
└────────────┬────────────────────┘
             │
             ▼
    CSV / XLSX / JSONL output
```

### Extraction strategy per field type

- **Bibliographic fields** (title, authors, year, journal, DOI): extracted from the first page using layout-aware scoring and noise filtering.
- **Organisms**: COMMON_SPECIES exact-match (weighted) + binomial regex, ranked by frequency.
- **Gene family**: regex patterns prioritised (e.g. `"WRKY gene family"`), fallback to known family name list.
- **Tools** (HMMER, BLAST, MEGA, etc.): exact keyword match with word-boundary anchors; BLAST additionally extracts the E-value.
- **Gene count**: multi-pattern regex covering "identified X genes", "a total of X", "comprises X members", etc.
- **Phylogenetic structure / duplication**: phrase-level matching for clades, subfamilies, tandem/segmental/WGD duplication.

---

## Extending the Tool

### Adding new gene families

Edit `KNOWN_GENE_FAMILIES` in `genomewide_pdf_extractor.py`:

```python
KNOWN_GENE_FAMILIES: list[str] = [
    ...,
    "YOUR_FAMILY",
]
```

### Adding new species

Edit `COMMON_SPECIES`:

```python
COMMON_SPECIES: list[str] = [
    ...,
    "Quercus robur",
]
```

### Adding new journals

Edit `JOURNAL_HINTS`:

```python
JOURNAL_HINTS: list[str] = [
    ...,
    "Your New Journal Name",
]
```

### Adding new output columns

1. Add the column name to `COLUMNS`.
2. Write an extractor method `extract_<fieldname>(self, text: str) -> str`.
3. Add a `_make_result(...)` call for your new field inside `extract_fields()`.

---

## Limitations

- Extraction accuracy depends on PDF text quality. Scanned/image-only PDFs will return mostly empty fields (consider running OCR first, e.g. with `ocrmypdf`).
- Multi-column PDF layouts may cause text to be extracted in reading order across columns, which can confuse section detection.
- The tool is tuned for plant genomics papers; animal/fungal studies are supported but the species and gene family lists are less comprehensive.

---

## Citation

If you use this tool in your research, please cite:

```
[Your Name] (2025). Genome-Wide PDF Extractor: A local tool for structured
metadata extraction from genome-wide gene family study PDFs.
GitHub: https://github.com/YOUR_USERNAME/genomewide-pdf-extractor
```

---

## License

MIT License — see [LICENSE](LICENSE) for full text.
