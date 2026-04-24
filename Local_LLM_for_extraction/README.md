# Local LLM Extraction Pipeline

This folder contains a standalone Python workflow for extracting structured metadata from research-paper PDFs related to genome-wide in silico studies.

## What the script does

`Complete_local_LLM.py` processes a directory of PDF papers and attempts to extract a standardized set of metadata columns, including:

- bibliographic fields such as title, authors, year, journal, DOI, and keywords
- biological context such as organism, gene family, and genome source
- methodological fields such as BLAST, HMMER, SMART, CDD, MEME, MCScanX, TBTools, qPCR, and RNA-seq usage
- summary findings such as number of genes identified, phylogenetic structure, and duplication events

The pipeline primarily uses deterministic heuristics over extracted PDF text. If a local Ollama model is available, it can optionally repair selected low-confidence fields.

## Core workflow

1. Extract text from PDFs with `pdfplumber` and/or `pypdf`
2. Score extraction quality and keep the stronger text output
3. Optionally pull table text with Camelot
4. Segment the paper into sections such as abstract, keywords, methods, results, and discussion
5. Apply rule-based extractors for each metadata field
6. Flag low-confidence values for manual review
7. Optionally ask a local Ollama model to improve selected weak fields

## Inputs

- A folder containing PDF research articles

## Outputs

For each run, the script writes:

- a main CSV with extracted metadata
- a review queue CSV for low-confidence cells
- a paper-level summary CSV
- an evidence JSONL file with value, confidence, and evidence snippets
- a run summary text file
- optionally an Excel copy of the main CSV

## Usage

Example:

```bash
python Complete_local_LLM.py --input path/to/pdfs --output-dir extractor_output --run-name local_extraction --excel
```

Optional arguments:

- `--review-threshold`: low-confidence cutoff, default `0.65`
- `--use-tables`: enables Camelot table extraction
- `--ollama-model`: local Ollama model name such as `qwen2.5:7b`

## Python dependencies

Required:

- `pandas`
- either `pdfplumber` or `pypdf`

Optional:

- `camelot` for table extraction
- a local Ollama server for low-confidence repair

## Privacy and deployment model

This pipeline is designed for local execution. PDF processing and optional LLM-assisted repair are intended to stay on the researcher's own machine. No remote API dependency is required by default.

## Project note

The extraction code is public, but the curated study dataset is not distributed in this repository. Data associated with the broader project are available upon request from the repository owner.
