import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".pydeps"))
os.environ["MPLCONFIGDIR"] = os.path.join(os.path.dirname(__file__), ".mplconfig")

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import textwrap

from repo_paths import configure_script_environment, resolve_data_file, resolve_output_dir

configure_script_environment()

DATA_FILE = resolve_data_file("Data22.csv", "Data_cleaned.csv", "Data.csv")
OUT_DIR = resolve_output_dir()
OUT_FILE = os.path.join(OUT_DIR, "Journal_Radial_StudyType_Figure.png")


JOURNAL_NORMALIZATION = {
    "PloS ONE": "PLOS ONE",
    "Plos One": "PLOS ONE",
}

TYPE_COLORS = {
    "In silico + RNA-seq + qPCR": "#cbd5e8",
    "In silico + RNA-seq": "#b3e2cd",
    "In silico + qPCR": "#fdcdac",
    "Only in silico": "#f4cae4",
}

RANK_COLORS = [
    "#a1c9f4", "#ffb482", "#8de5a1", "#ff9f9b", "#d0bbff", "#debb9b",
    "#fab0e4", "#cfcfcf", "#fffea3", "#b9f2f0", "#e8e5a1", "#f2c6b3",
]


def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)


def clean_text(value):
    text = "" if pd.isna(value) else str(value).strip()
    text = " ".join(text.split())
    return JOURNAL_NORMALIZATION.get(text, text)


def yes(value):
    text = "" if pd.isna(value) else str(value).strip().lower()
    return text not in {"", "no", "none", "false", "n/a", "na"} and "not used" not in text


def compact(text, n=30):
    return text if len(text) <= n else text[: n - 3] + "..."


def classify_study(row):
    q = yes(row["Experimental Validation (qPCR)"])
    r = yes(row["Expression Analysis (RNA-seq)"])
    if q and r:
        return "In silico + RNA-seq + qPCR"
    if q:
        return "In silico + qPCR"
    if r:
        return "In silico + RNA-seq"
    return "Only in silico"


def main():
    ensure_dir(OUT_DIR)

    df = pd.read_csv(DATA_FILE, encoding="latin-1")
    df["Journal_clean"] = df["Journal"].map(clean_text).replace("", "Unknown")
    df["Study_Type"] = df.apply(classify_study, axis=1)

    top_journals = df["Journal_clean"].value_counts().head(12).index.tolist()
    counts = df["Journal_clean"].value_counts().reindex(top_journals)
    comp = (
        df[df["Journal_clean"].isin(top_journals)]
        .groupby(["Journal_clean", "Study_Type"])
        .size()
        .unstack(fill_value=0)
        .reindex(top_journals)
        .reindex(columns=list(TYPE_COLORS.keys()), fill_value=0)
    )
    comp_pct = comp.div(comp.sum(axis=1), axis=0)

    n = len(top_journals)
    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
    width = 2 * np.pi / n * 0.9

    fig = plt.figure(figsize=(14, 12), facecolor="white")
    ax = plt.subplot(111, polar=True)
    ax.set_facecolor("white")
    ax.set_theta_direction(-1)
    ax.set_theta_offset(np.pi / 2)
    ax.set_ylim(0, 1.25)
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["polar"].set_visible(False)

    max_count = counts.max()
    inner_base = 0.38
    inner_thickness = 0.18
    outer_base = 0.60
    outer_max = 0.34

    # subtle ring guides
    guide_t = np.linspace(0, 2 * np.pi, 400)
    ax.plot(guide_t, np.full_like(guide_t, inner_base), color="#e7eef6", linewidth=1)
    ax.plot(guide_t, np.full_like(guide_t, outer_base), color="#edf3f8", linewidth=1)

    for i, journal in enumerate(top_journals):
        angle = theta[i]

        # Inner stacked study-type ring
        start = inner_base
        for category, color in TYPE_COLORS.items():
            h = inner_thickness * comp_pct.loc[journal, category]
            if h > 0:
                ax.bar(angle, h, width=width * 0.86, bottom=start, color=color, edgecolor="white", linewidth=1.0, align="center", zorder=3)
                start += h

        # Outer ranking bar
        bar_h = outer_max * counts.loc[journal] / max_count
        ax.bar(
            angle,
            bar_h,
            width=width * 0.92,
            bottom=outer_base,
            color=RANK_COLORS[i % len(RANK_COLORS)],
            edgecolor="white",
            linewidth=1.2,
            align="center",
            zorder=4,
        )

        # Count label on outer bar
        ax.text(
            angle,
            outer_base + bar_h + 0.035,
            str(int(counts.loc[journal])),
            ha="center",
            va="center",
            fontsize=9,
            color="#16324f",
            fontweight="bold",
            rotation=np.degrees(-angle + np.pi / 2),
            rotation_mode="anchor",
        )

    # Radially placed wrapped labels directly outside numbers
    for i, journal in enumerate(top_journals):
        angle = theta[i]
        bar_h = outer_max * counts.loc[journal] / max_count
        
        # Calculate screen angle for proper text orientation
        screen_angle_deg = (90 - np.degrees(angle)) % 360
        
        if 90 < screen_angle_deg <= 270:
            rot = screen_angle_deg + 180
            ha = "right"
        else:
            rot = screen_angle_deg
            ha = "left"
            
        journal_wrapped = textwrap.fill(journal, width=18)
        
        ax.text(
            angle,
            outer_base + bar_h + 0.07,
            journal_wrapped,
            ha=ha,
            va="center",
            fontsize=12.5,
            fontweight="bold",
            color="#1f2937",
            rotation=rot,
            rotation_mode="anchor",
            zorder=6
        )

    # center annotation
    ax.text(0, 0.08, "Top publishing\njournals in the\ndataset", ha="center", va="center", fontsize=15, fontweight="bold", color="#0b2540")

    # legend
    legend_handles = [
        plt.Line2D([0], [0], color=color, lw=8)
        for color in TYPE_COLORS.values()
    ]
    fig.legend(
        legend_handles,
        list(TYPE_COLORS.keys()),
        loc="lower center",
        ncol=2,
        frameon=False,
        fontsize=10,
        bbox_to_anchor=(0.5, 0.04),
    )

    plt.savefig(OUT_FILE, dpi=320, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print("Journal radial study-type figure generated.")
    print("Output:", OUT_FILE)


if __name__ == "__main__":
    main()
