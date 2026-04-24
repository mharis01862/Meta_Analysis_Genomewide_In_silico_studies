import os
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
DATA_DIR = PROJECT_ROOT / "Data"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "Generated_Figures"


def configure_script_environment() -> None:
    sys.path.insert(0, str(SCRIPT_DIR / ".pydeps"))
    os.environ.setdefault("MPLCONFIGDIR", str(SCRIPT_DIR / ".mplconfig"))


def resolve_data_file(*candidates: str) -> str:
    override = os.getenv("META_ANALYSIS_DATA_FILE")
    if override:
        return override

    search_roots = [PROJECT_ROOT, DATA_DIR, SCRIPT_DIR]
    for candidate in candidates:
        candidate_path = Path(candidate)
        if candidate_path.is_absolute() and candidate_path.exists():
            return str(candidate_path)
        for root in search_roots:
            path = root / candidate
            if path.exists():
                return str(path)

    fallback = DATA_DIR / (candidates[0] if candidates else "Data.csv")
    return str(fallback)


def resolve_output_dir(*parts: str) -> str:
    base = Path(os.getenv("META_ANALYSIS_OUTPUT_DIR", str(DEFAULT_OUTPUT_DIR)))
    return str(base.joinpath(*parts))


def resolve_auxiliary_file(filename: str) -> str:
    for root in (SCRIPT_DIR, PROJECT_ROOT, DATA_DIR):
        path = root / filename
        if path.exists():
            return str(path)
    return str(PROJECT_ROOT / filename)
