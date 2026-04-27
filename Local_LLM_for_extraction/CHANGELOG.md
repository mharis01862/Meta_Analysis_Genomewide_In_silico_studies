# Changelog

All notable changes are documented here following [Keep a Changelog](https://keepachangelog.com/) conventions.

---

## [1.1.0] — 2025

### Added
- `--file` flag: process a single PDF instead of a whole directory
- `--verbose` flag: print every extracted field + confidence score to stdout
- `process_single_file()` method on `GenomicsExtractor`
- WGD (whole-genome duplication) detection in `extract_duplication()`
- 7 gene-count regex patterns (was 4): added "comprises", "contains", "consists of"
- 13 new `KNOWN_GENE_FAMILIES` entries (kinases, transporters, stress-response families)
- 15 new `COMMON_SPECIES` entries (cotton, banana, quinoa, spinach, and others)
- 7 new `JOURNAL_HINTS` entries (Plant Cell, Planta, Crop Science, and others)
- Additional genome databases: RAP-DB, CoGe, Sol Genomics, CottonFGD, PLAZA
- LLM repair now covers 9 fields (was 7): added Country of Origin and Year
- `--verbose` detailed per-field output for debugging single papers

### Changed
- Class renamed `CompleteLocalExtractor` → `GenomicsExtractor` (cleaner public API)
- All private helpers prefixed with `_` to separate public from internal API
- `_reader_score()` adds +10 bonus for "genome-wide" signal on first page
- Title extraction scores +8 for "transcription factor" / "protein kinase" patterns
- BLAST E-value regex handles `×10⁻⁵` Unicode multiplication notation
- Author-line cleaning extended to strip `‡ § ¶` symbols in addition to `* # †`
- Confidence scoring improved for HMMER (Pfam accession present → 0.95),
  Gene Family (known family match → +0.15), Organism (known species → +0.15)
- `process_directory()` prints organism | gene family | gene count per paper as QC
- Full module-level docstring with embedded MIT licence text

### Fixed
- f-string with `+` concatenation inside conditional (latent SyntaxWarning on Python 3.12)
- `_extract_other_tools` now uses word-boundary regex instead of `str.lower() in`
  to avoid false positives (e.g. "smart" matching inside "smartwatch")

---

## [1.0.0] — 2025

### Initial release
- Dual PDF reader back-end (pdfplumber + pypdf) with automatic quality selection
- 30-field structured extraction for genome-wide gene family identification papers
- Section-aware extraction (Abstract, Methods, Results, Discussion)
- Confidence scoring (0–1) per field with review-queue output
- Evidence JSONL for full field-level provenance
- Optional Camelot table extraction (`--use-tables`)
- Optional Ollama local LLM repair (`--ollama-model`)
- Excel output (`--excel`)
- Batch directory processing with recursive PDF discovery
