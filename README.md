# Meta-Analysis of Genome-Wide In Silico Studies

This repository contains the code, generated figures, and supporting project materials for a meta-analysis focused on genome-wide in silico studies. The public repository is intentionally limited to publishable project assets only.

The analysis dataset and manuscript source files are not public in this repository. The underlying data contain non-public research materials and are available upon reasonable request from the project owner.

## What is included

- Figure-generation scripts used to produce the visual summaries in the study
- Rendered figure files and figure source documents
- A local large-language-model extraction pipeline for structured paper metadata extraction
- Project-level documentation for reproducing the code-side workflow once data access is granted

## What is not included

- `Data/`
- `Manuscript/`

Both directories are excluded from version control for privacy and publication-control reasons.

## Repository structure

- `Figure_Generation_Codes/`: Python scripts used to generate the study figures
- `Figures/`: exported figure panels used in the project
- `Figure_seperate_files/`: individual figure documents
- `Merged_figures/`: merged figure composites
- `Local_LLM_for_extraction/`: standalone local PDF extraction pipeline and its documentation

## Data access policy

The study data are not publicly distributed through this repository. Data access is available upon request. Please contact the repository owner to request the curated dataset used for the analysis.

When data are provided, place the supplied CSV file(s) inside `Data/` and review the script-specific notes below.

## Figure-generation code

The scripts in `Figure_Generation_Codes/` were adjusted for repository use so they no longer depend on the original author machine paths. They now resolve data and output locations relative to this repository and can also be configured with environment variables.

### Default behavior

- Input data are searched in `Data/`
- Outputs are written to `Generated_Figures/`
- Matplotlib cache data are written inside each script folder's local `.mplconfig/`

### Optional environment variables

- `META_ANALYSIS_DATA_FILE`: explicit path to a CSV file to use
- `META_ANALYSIS_OUTPUT_DIR`: explicit directory for figure exports

### Data-file expectations

Some scripts historically referenced `Data.csv`, `Data_cleaned.csv`, or `Data22.csv`. The repository helper now tries multiple common names, but the public repo does not include those files. If you receive the private dataset, keep the original column names intact.

### Known external dependency

`Generate_Country_Origin_Map.py` expects a `world_countries.geojson` file for map rendering. That auxiliary file is not present in this public snapshot and should be supplied separately when reproducing the map.

## Local LLM extraction workflow

The extraction pipeline in `Local_LLM_for_extraction/` is documented separately in `Local_LLM_for_extraction/README.md`. In short, it:

- extracts text from research-paper PDFs
- applies rule-based heuristics to recover structured metadata
- optionally uses a local Ollama model to repair low-confidence fields
- emits CSV, review-queue, summary, and evidence files for manual validation

This component is also released under an MIT license for code only.

## Reproducibility notes

This repository is intended as a code and figure companion to the study, not as a public release of the private working dataset. Full regeneration of all outputs requires access to the non-public data files and, for some scripts, auxiliary assets that are not bundled here.

## License

The code in this repository is released under the MIT License. The MIT license applies to code only. It does not grant rights to unpublished data, manuscript drafts, or any third-party content that may be referenced by the project.
