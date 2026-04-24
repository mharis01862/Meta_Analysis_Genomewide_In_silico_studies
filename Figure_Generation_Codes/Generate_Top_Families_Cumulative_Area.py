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
from matplotlib.patches import Rectangle

from repo_paths import configure_script_environment, resolve_data_file, resolve_output_dir

configure_script_environment()

DATA_FILE = resolve_data_file("Data22.csv", "Data_cleaned.csv", "Data.csv")
OUT_DIR = resolve_output_dir("Figure_07_Top_Families_Cumulative_Area")
OUT_FILE = os.path.join(OUT_DIR, "Top_9_Families_Cumulative_Area.png")


FAMILY_COLORS = [
    "#45275a",
    "#5b3bb3",
    "#13d7ba",
    "#f0a23d",
    "#f7c34a",
    "#c95f98",
    "#8d1209",
    "#0b4f8c",
    "#6a4c93",
]


def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)


def clean_text(value):
    text = "" if pd.isna(value) else str(value).strip()
    return " ".join(text.split())


def parse_year(value):
    m = re.search(r"(19|20)\d{2}", clean_text(value))
    return int(m.group(0)) if m else np.nan


def normalize_family(value):
    text = clean_text(value)
    low = text.lower()
    if "wrky" in low:
        return "WRKY"
    if "dof" in low:
        return "Dof"
    if "hsf" in low or "heat shock" in low:
        return "HSF"
    if "squamosa" in low or re.fullmatch(r"spl.*", low):
        return "SPL"
    if "sweet" in low:
        return "SWEET"
    if "mapkkk" in low or "mapk" in low:
        return "MAPKKK"
    if "bzip" in low:
        return "bZIP"
    if "sod" in low or "superoxide" in low:
        return "SOD"
    if "ap2" in low or "erf" in low or "erebp" in low:
        return "AP2/ERF"
    if "gst" in low:
        return "GST"
    if "gata" in low:
        return "GATA"
    if "expansin" in low:
        return "Expansin"
    if "aquaporin" in low:
        return "Aquaporin"
    if "f-box" in low or "f box" in low:
        return "F-box"
    if "wox" in low:
        return "WOX"
    if "hsp20" in low:
        return "Hsp20"
    if "terpene synthase" in low or low == "tps":
        return "TPS"
    return text


def main():
    ensure_dir(OUT_DIR)

    df = pd.read_csv(DATA_FILE, encoding="latin-1")
    df["Year_num"] = df["Year"].map(parse_year)
    df["Family_norm"] = df["Gene Family"].map(normalize_family)
    df = df.dropna(subset=["Year_num"])
    df["Year_num"] = df["Year_num"].astype(int)

    family_order = (
        df["Family_norm"]
        .value_counts()
        .head(8)
        .index
        .tolist()
    )

    all_years = np.arange(df["Year_num"].min(), df["Year_num"].max() + 1)
    family_year = (
        df[df["Family_norm"].isin(family_order)]
        .groupby(["Family_norm", "Year_num"])
        .size()
        .rename("n")
        .reset_index()
    )

    fig, axes = plt.subplots(2, 4, figsize=(15, 7.9), facecolor="white", sharex=True)
    axes = axes.flatten()

    max_cum = 1
    family_series = {}
    for family in family_order:
        sub = family_year[family_year["Family_norm"] == family][["Year_num", "n"]].set_index("Year_num").reindex(all_years, fill_value=0)
        sub["cum"] = sub["n"].cumsum()
        family_series[family] = sub
        max_cum = max(max_cum, int(sub["cum"].max()))

    for idx, family in enumerate(family_order):
        ax = axes[idx]
        ax.set_facecolor("white")
        color = FAMILY_COLORS[idx % len(FAMILY_COLORS)]
        sub = family_series[family]

        ax.fill_between(all_years, sub["cum"].values, color=color, alpha=0.78, linewidth=0)
        ax.plot(all_years, sub["cum"].values, color=color, linewidth=2.2)

        ax.grid(axis="both", color="#dfe7f1", linewidth=1, alpha=0.7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.tick_params(axis="both", length=0, labelsize=10, colors="#5a6d7f")
        ax.set_ylim(-0.2, max_cum + 0.2)
        ax.set_yticks(np.arange(0, max_cum + 1, 1 if max_cum <= 8 else 2))

        # Header strip per panel
        ax.add_patch(
            Rectangle(
                (0, 1.02),
                1,
                0.12,
                transform=ax.transAxes,
                facecolor="#3f556c",
                edgecolor="#2f4255",
                linewidth=1.2,
                clip_on=False,
                zorder=10,
            )
        )
        ax.text(
            0.5,
            1.08,
            family,
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=12.5,
            color="white",
            fontweight="bold",
            zorder=11,
        )

    # Hide any unused axes if fewer than 8 families
    for j in range(len(family_order), len(axes)):
        axes[j].axis("off")

    for ax in axes[-4:]:
        ax.set_xticks(np.arange(all_years.min(), all_years.max() + 1, 5))
        ax.set_xticklabels([str(y) for y in np.arange(all_years.min(), all_years.max() + 1, 5)], fontsize=10.5)

    fig.suptitle(
        "Faceted Cumulative Area Chart of Top Gene Families",
        x=0.07,
        y=0.975,
        ha="left",
        fontsize=21,
        fontweight="bold",
        color="#334a62",
    )
    fig.text(0.03, 0.50, "Total Publications (Cumulative)", rotation=90, va="center", ha="center", fontsize=16, color="#1f2937")
    fig.text(0.50, 0.055, "Year", va="center", ha="center", fontsize=17, color="#1f2937")

    plt.tight_layout(rect=[0.05, 0.08, 0.995, 0.95], w_pad=0.6, h_pad=1.0)
    plt.savefig(OUT_FILE, dpi=500, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print("Top families cumulative area figure generated.")
    print("Output:", OUT_FILE)


if __name__ == "__main__":
    main()
