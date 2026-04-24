import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".pydeps"))
os.environ["MPLCONFIGDIR"] = os.path.join(os.path.dirname(__file__), ".mplconfig")

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch

from repo_paths import configure_script_environment, resolve_data_file, resolve_output_dir

configure_script_environment()

DATA_FILE = resolve_data_file("Data22.csv", "Data_cleaned.csv", "Data.csv")
OUT_DIR = resolve_output_dir("Figure_06_Top_Families_Radar")
OUT_FILE = os.path.join(OUT_DIR, "Top_Families_Capability_Radar.png")


AXIS_LABELS = ["Domain", "Motifs", "Phylogeny", "Synteny", "Validation"]
FAMILY_ORDER = ["WRKY", "Dof", "HSF", "SPL", "SWEET", "MAPKKK", "bZIP", "SOD"]
FAMILY_COLORS = {
    "WRKY": "#43215e",
    "Dof": "#0fd6b2",
    "HSF": "#f5b21d",
    "SPL": "#8f120b",
    "SWEET": "#0b4f8c",
    "MAPKKK": "#e76f51",
    "bZIP": "#6a4c93",
    "SOD": "#2a9d8f",
}


def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)


def clean_text(value):
    text = "" if pd.isna(value) else str(value).strip()
    return " ".join(text.split())


def yes(value):
    text = clean_text(value).lower()
    return text not in {"", "no", "none", "false", "n/a", "na"} and "not used" not in text


def normalize_family(value):
    text = clean_text(value)
    low = text.lower()
    if "wrky" in low:
        return "WRKY"
    if "dof" in low:
        return "Dof"
    if "squamosa" in low or re.fullmatch(r"spl.*", low):
        return "SPL"
    if "hsf" in low or "heat shock" in low:
        return "HSF"
    if "sweet" in low:
        return "SWEET"
    if "mapkkk" in low or "mapk" in low:
        return "MAPKKK"
    if "bzip" in low:
        return "bZIP"
    if "sod" in low or "superoxide" in low:
        return "SOD"
    return text


def family_scores(df):
    df["Family_norm"] = df["Gene Family"].map(normalize_family)
    df["Capability_Domain"] = df[[
        "HMMER (with Pfam domains)",
        "BLAST (E-value cut-off used)",
        "SMART",
        "CDD",
    ]].apply(lambda row: any(yes(v) for v in row), axis=1)
    df["Capability_Motifs"] = df[[
        "MEME",
        "Others (Tools for Motif Analysis)",
    ]].apply(lambda row: any(yes(v) for v in row), axis=1)
    df["Capability_Phylogeny"] = df[[
        "MEGA",
        "IQ-TREE",
        "Phylogenetic Structure",
    ]].apply(lambda row: any(yes(v) for v in row), axis=1)
    df["Capability_Synteny"] = df[[
        "MCScanX",
        "Other (Tools for Synteny)",
    ]].apply(lambda row: any(yes(v) for v in row), axis=1)
    df["Capability_Validation"] = df[[
        "Experimental Validation (qPCR)",
        "Expression Analysis (RNA-seq)",
    ]].apply(lambda row: any(yes(v) for v in row), axis=1)

    scores = (
        df[df["Family_norm"].isin(FAMILY_ORDER)]
        .groupby("Family_norm")[[
            "Capability_Domain",
            "Capability_Motifs",
            "Capability_Phylogeny",
            "Capability_Synteny",
            "Capability_Validation",
        ]]
        .mean()
        .reindex(FAMILY_ORDER)
        .fillna(0)
    )
    counts = df["Family_norm"].value_counts().reindex(FAMILY_ORDER).fillna(0).astype(int)
    return scores, counts


def add_family_header(ax, family):
    box = FancyBboxPatch(
        (-0.06, 1.06),
        1.12,
        0.12,
        transform=ax.transAxes,
        boxstyle="round,pad=0.02,rounding_size=0.02",
        facecolor="#3f556c",
        edgecolor="#31465a",
        linewidth=1.0,
        clip_on=False,
        zorder=10,
    )
    ax.add_patch(box)
    ax.text(
        0.5,
        1.12,
        family,
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=11.2,
        color="white",
        fontweight="bold",
        zorder=11,
    )
def main():
    ensure_dir(OUT_DIR)
    df = pd.read_csv(DATA_FILE, encoding="latin-1")
    scores, counts = family_scores(df)

    angles = np.linspace(0, 2 * np.pi, len(AXIS_LABELS), endpoint=False).tolist()
    angles += angles[:1]

    fig, axes = plt.subplots(2, 4, figsize=(14.5, 8.6), subplot_kw=dict(polar=True), facecolor="white")
    axes = axes.flatten()

    for idx, (ax, family) in enumerate(zip(axes, FAMILY_ORDER)):
        vals = scores.loc[family].tolist()
        vals += vals[:1]
        color = FAMILY_COLORS[family]

        ax.set_facecolor("white")
        ax.set_theta_offset(np.pi / 2)
        ax.set_theta_direction(-1)

        # Glow and main profile
        ax.plot(angles, vals, color=color, linewidth=5.4, alpha=0.08)
        ax.plot(angles, vals, color=color, linewidth=2.1)
        ax.fill(angles, vals, color=color, alpha=0.18)

        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(AXIS_LABELS, fontsize=7.6, color="#5f7182")
        ax.set_yticks([0.25, 0.5, 0.75, 1.0])
        ax.set_yticklabels([])
        ax.set_ylim(0, 1.0)
        ax.grid(color="#b7cadf", linestyle=(0, (3, 4)), linewidth=1.15)
        ax.spines["polar"].set_visible(False)
        add_family_header(ax, family)
        count_y = -0.06 if idx < 4 else -0.11
        ax.text(
            0.5,
            count_y,
            f"n={counts.loc[family]}",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=9.6,
            color="#3f556c",
            fontweight="bold",
        )

    fig.suptitle(
        "Capability Radar by Top Gene Families",
        x=0.07,
        y=0.975,
        ha="left",
        fontsize=21,
        fontweight="bold",
        color="#27415d",
    )
    fig.text(0.03, 0.50, "Score", rotation=90, va="center", ha="center", fontsize=13, color="#24384d")
    fig.text(0.50, 0.055, "Capability", va="center", ha="center", fontsize=14, color="#24384d")

    plt.tight_layout(rect=[0.04, 0.08, 0.99, 0.95], w_pad=1.0, h_pad=1.0)
    plt.savefig(OUT_FILE, dpi=500, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print("Top families radar figure generated.")
    print("Output:", OUT_FILE)


if __name__ == "__main__":
    main()
