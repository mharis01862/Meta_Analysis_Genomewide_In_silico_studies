import csv
import os
import re
from collections import Counter, defaultdict

import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".pydeps"))
os.environ["MPLCONFIGDIR"] = os.path.join(os.path.dirname(__file__), ".mplconfig")

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Rectangle

from repo_paths import configure_script_environment, resolve_data_file, resolve_output_dir

configure_script_environment()

DATA_FILE = resolve_data_file("Data_cleaned.csv", "Data22.csv", "Data.csv")
OUT_DIR = resolve_output_dir("Figure_08_Genome_Source_Ecosystem")
OUT_FILE = os.path.join(OUT_DIR, "Genome_Source_Ecosystem_Flow.png")


SOURCE_COLORS = {
    "Phytozome": "#0b4f8c",
    "NCBI": "#2a9d8f",
    "Ensembl Plants": "#5b3bb3",
    "Sol Genomics Network": "#e76f51",
    "CottonGen": "#f4a261",
    "TAIR": "#6a4c93",
    "GDR": "#c95f98",
    "Other sources": "#9aa9b8",
}

TYPE_COLORS = {
    "In silico + RNA-seq + qPCR": "#0b4f8c",
    "In silico + RNA-seq": "#2a9d8f",
    "In silico + qPCR": "#f4a261",
    "Only in silico": "#d9e6f2",
}


def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)


def clean_text(text):
    text = "" if text is None else str(text).strip()
    text = re.sub(r"\s+", " ", text)
    return text.strip(" ,;")


def compact(text, n=24):
    text = clean_text(text)
    return text if len(text) <= n else text[: n - 3] + "..."


def yes(value):
    t = clean_text(value).lower()
    return t not in {"", "no", "none", "false", "n/a", "na"} and "not used" not in t


def normalize_source(value):
    text = clean_text(value)
    low = text.lower()
    if not text:
        return "Unspecified"
    if "phytozome" in low:
        return "Phytozome"
    if "ensembl" in low:
        return "Ensembl Plants"
    if "ncbi" in low:
        return "NCBI"
    if "cotton" in low:
        return "CottonGen"
    if "sol genomics" in low or "sgn" in low:
        return "Sol Genomics Network"
    if "gramene" in low:
        return "Gramene"
    if "tair" in low:
        return "TAIR"
    if "gdr" in low:
        return "GDR"
    if "national genomics data center" in low:
        return "NGDC"
    if "joint genome institute" in low or "jgi" in low:
        return "JGI"
    return text


def classify_type(row):
    q = yes(row.get("Experimental Validation (qPCR)", ""))
    r = yes(row.get("Expression Analysis (RNA-seq)", ""))
    if q and r:
        return "In silico + RNA-seq + qPCR"
    if q:
        return "In silico + qPCR"
    if r:
        return "In silico + RNA-seq"
    return "Only in silico"


def load_rows():
    with open(DATA_FILE, encoding="latin-1", newline="") as f:
        rows = list(csv.DictReader(f))
    normalized = []
    source_counts = Counter()
    organism_counts = Counter()
    for row in rows:
        source = normalize_source(row.get("Genome Source", ""))
        organism = clean_text(row.get("Primary Organism", "")) or "Unspecified"
        study_type = classify_type(row)
        normalized.append((source, organism, study_type))
        source_counts[source] += 1
        organism_counts[organism] += 1
    return normalized, source_counts, organism_counts


def stack_positions(nodes, counts, total_height=0.84, start_y=0.92, gap=0.014):
    total = sum(counts[n] for n in nodes)
    usable = total_height - gap * (len(nodes) - 1)
    positions = {}
    cursor = start_y
    for node in nodes:
        h = usable * counts[node] / total if total else 0
        positions[node] = [cursor - h, cursor]
        cursor = cursor - h - gap
    return positions


def add_flow(ax, x0, x1, y0a, y0b, y1a, y1b, color, alpha=0.5):
    ctrl = 0.16
    verts = [
        (x0, y0a),
        (x0 + ctrl, y0a),
        (x1 - ctrl, y1a),
        (x1, y1a),
        (x1, y1b),
        (x1 - ctrl, y1b),
        (x0 + ctrl, y0b),
        (x0, y0b),
        (x0, y0a),
    ]
    codes = [
        Path.MOVETO,
        Path.CURVE4,
        Path.CURVE4,
        Path.CURVE4,
        Path.LINETO,
        Path.CURVE4,
        Path.CURVE4,
        Path.CURVE4,
        Path.CLOSEPOLY,
    ]
    ax.add_patch(PathPatch(Path(verts, codes), facecolor=color, edgecolor="none", alpha=alpha, zorder=1))


def main():
    ensure_dir(OUT_DIR)
    rows, source_counts_all, organism_counts_all = load_rows()

    top_sources = [s for s, _ in source_counts_all.most_common(7)]
    top_organisms = [o for o, _ in organism_counts_all.most_common(8)]

    records = []
    source_counts = Counter()
    organism_counts = Counter()
    type_counts = Counter()
    for source, organism, study_type in rows:
        source_group = source if source in top_sources else "Other sources"
        organism_group = organism if organism in top_organisms else "Other organisms"
        records.append((source_group, organism_group, study_type))
        source_counts[source_group] += 1
        organism_counts[organism_group] += 1
        type_counts[study_type] += 1

    sources = [s for s in top_sources if source_counts[s] > 0]
    if source_counts["Other sources"] > 0:
        sources.append("Other sources")
    organisms = [o for o in top_organisms if organism_counts[o] > 0]
    if organism_counts["Other organisms"] > 0:
        organisms.append("Other organisms")
    types = ["Only in silico", "In silico + qPCR", "In silico + RNA-seq", "In silico + RNA-seq + qPCR"]

    source_pos = stack_positions(sources, source_counts)
    organism_pos = stack_positions(organisms, organism_counts)
    type_pos = stack_positions(types, type_counts)

    src_org = Counter((s, o) for s, o, _ in records)
    org_type = Counter((o, t) for _, o, t in records)

    src_offsets = {s: source_pos[s][0] for s in sources}
    org_in_offsets = {o: organism_pos[o][0] for o in organisms}
    org_out_offsets = {o: organism_pos[o][0] for o in organisms}
    type_offsets = {t: type_pos[t][0] for t in types}

    fig, ax = plt.subplots(figsize=(16.5, 9.5), facecolor="white")
    ax.set_facecolor("white")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    x_source, x_org, x_type = 0.12, 0.50, 0.88
    bar_w = 0.028

    # flows source -> organism
    for source in sources:
        for organism in organisms:
            n = src_org.get((source, organism), 0)
            if n <= 0:
                continue
            y0a = src_offsets[source]
            y0b = y0a + (source_pos[source][1] - source_pos[source][0]) * n / source_counts[source]
            y1a = org_in_offsets[organism]
            y1b = y1a + (organism_pos[organism][1] - organism_pos[organism][0]) * n / organism_counts[organism]
            add_flow(ax, x_source + bar_w, x_org - bar_w, y0a, y0b, y1a, y1b, SOURCE_COLORS.get(source, "#9aa9b8"), alpha=0.34)
            src_offsets[source] = y0b
            org_in_offsets[organism] = y1b

    # flows organism -> study type
    for organism in organisms:
        for study_type in types:
            n = org_type.get((organism, study_type), 0)
            if n <= 0:
                continue
            y0a = org_out_offsets[organism]
            y0b = y0a + (organism_pos[organism][1] - organism_pos[organism][0]) * n / organism_counts[organism]
            y1a = type_offsets[study_type]
            y1b = y1a + (type_pos[study_type][1] - type_pos[study_type][0]) * n / type_counts[study_type]
            add_flow(ax, x_org + bar_w, x_type - bar_w, y0a, y0b, y1a, y1b, TYPE_COLORS[study_type], alpha=0.36)
            org_out_offsets[organism] = y0b
            type_offsets[study_type] = y1b

    # node bars and labels
    for source in sources:
        y0, y1 = source_pos[source]
        ax.add_patch(Rectangle((x_source, y0), bar_w, y1 - y0, facecolor=SOURCE_COLORS.get(source, "#9aa9b8"), edgecolor="white", linewidth=1.2, zorder=3))
        cy = (y0 + y1) / 2
        source_label_offsets = {
            "CottonGen": (0.012, -0.009),
            "TAIR": (0.009, -0.006),
            "GDR": (0.006, -0.008),
        }
        name_dy, count_dy = source_label_offsets.get(source, (0.010, -0.012))
        ax.text(x_source - 0.015, cy + name_dy, compact(source, 22), ha="right", va="center", fontsize=10.5, color="#16324f", fontweight="bold")
        ax.text(x_source - 0.015, cy + count_dy, f"n={source_counts[source]}", ha="right", va="center", fontsize=8.5, color="#5a6d7f")

    for organism in organisms:
        y0, y1 = organism_pos[organism]
        ax.add_patch(Rectangle((x_org - bar_w / 2, y0), bar_w, y1 - y0, facecolor="#d9e6f2", edgecolor="#c8d8e7", linewidth=1.0, zorder=3))
        # Italicize scientific names (exclude generic labels)
        is_scientific = organism not in ["Other organisms", "Unspecified"]
        font_style = "italic" if is_scientific else "normal"
        ax.text(x_org, y1 + 0.008, compact(organism, 24), ha="center", va="bottom", fontsize=9.2, color="#16324f", fontweight="bold", fontstyle=font_style)

    for study_type in types:
        y0, y1 = type_pos[study_type]
        ax.add_patch(Rectangle((x_type - bar_w, y0), bar_w, y1 - y0, facecolor=TYPE_COLORS[study_type], edgecolor="white", linewidth=1.2, zorder=3))
        cy = (y0 + y1) / 2
        ax.text(x_type + 0.015, cy + 0.010, study_type, ha="left", va="center", fontsize=10.2, color="#16324f", fontweight="bold")
        ax.text(x_type + 0.015, cy - 0.012, f"n={type_counts[study_type]}", ha="left", va="center", fontsize=8.5, color="#5a6d7f")

    ax.text(x_source + bar_w / 2, 0.965, "Genome Sources", ha="center", va="center", fontsize=14, color="#0b4f8c", fontweight="bold")
    ax.text(x_org, 0.965, "Primary Organisms", ha="center", va="center", fontsize=14, color="#4a6179", fontweight="bold")
    ax.text(x_type - bar_w / 2, 0.965, "Study Design", ha="center", va="center", fontsize=14, color="#2a9d8f", fontweight="bold")

    ax.set_title(
        "Genome Source Ecosystem Flow",
        loc="left",
        fontsize=22,
        pad=18,
        color="#0b2540",
        fontweight="bold",
    )

    plt.tight_layout()
    plt.savefig(OUT_FILE, dpi=450, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print("Genome source ecosystem figure generated.")
    print("Output:", OUT_FILE)


if __name__ == "__main__":
    main()
