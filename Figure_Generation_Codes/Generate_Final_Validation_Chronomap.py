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

from repo_paths import configure_script_environment, resolve_data_file, resolve_output_dir

configure_script_environment()

DATA_FILE = resolve_data_file("Data.csv", "Data_cleaned.csv", "Data22.csv")
OUT_DIR = resolve_output_dir("Figure_13_Validation_Expression")
OUT_FILE = os.path.join(OUT_DIR, "Validation_Expression_Chronomap.png")

STATE_ORDER = [
    "Only in silico",
    "qPCR only",
    "RNA-seq only",
    "qPCR + RNA-seq",
]

STATE_COLORS = {
    "Only in silico": "#d7e8f6",
    "qPCR only": "#f2c14e",
    "RNA-seq only": "#58b89c",
    "qPCR + RNA-seq": "#1f7a8c",
}


def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)


def clean_text(value):
    if pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def yes(value):
    text = clean_text(value).lower()
    return text not in {"", "no", "none", "false", "n/a", "na"} and "not used" not in text


def parse_year(value):
    match = re.search(r"(19|20)\d{2}", clean_text(value))
    return int(match.group(0)) if match else np.nan


def classify_state(row):
    qpcr = yes(row["Experimental Validation (qPCR)"])
    rnaseq = yes(row["Expression Analysis (RNA-seq)"])
    if qpcr and rnaseq:
        return "qPCR + RNA-seq"
    if qpcr:
        return "qPCR only"
    if rnaseq:
        return "RNA-seq only"
    return "Only in silico"


def main():
    ensure_dir(OUT_DIR)

    df = pd.read_csv(DATA_FILE, encoding="latin-1")
    df["Year_num"] = df["Year"].map(parse_year)
    df = df.dropna(subset=["Year_num"]).copy()
    df["Evidence_State"] = df.apply(classify_state, axis=1)

    years = list(range(int(df["Year_num"].min()), int(df["Year_num"].max()) + 1))
    counts = (
        df.groupby(["Year_num", "Evidence_State"]).size().reset_index(name="count")
    )

    angle_lookup = {
        year: angle
        for year, angle in zip(years, np.linspace(np.pi / 2, np.pi / 2 - 2 * np.pi, len(years), endpoint=False))
    }
    radius_lookup = {
        "Only in silico": 1.55,
        "qPCR only": 2.45,
        "RNA-seq only": 3.35,
        "qPCR + RNA-seq": 4.35,
    }

    fig = plt.figure(figsize=(14.8, 12.6), facecolor="#f4f7fb")
    ax = fig.add_subplot(111, projection="polar")
    ax.set_facecolor("#fbfdff")
    ax.set_theta_direction(-1)
    ax.set_theta_offset(np.pi / 2.0)
    ax.set_ylim(0, 5.15)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    ax.spines["polar"].set_visible(False)

    # Ring guides
    theta = np.linspace(0, 2 * np.pi, 500)
    for state, radius in radius_lookup.items():
        ax.plot(theta, np.full_like(theta, radius), color="#dfe8f1", linewidth=1.2, zorder=1)

    # Year spokes and labels
    for year in years:
        ang = angle_lookup[year]
        ax.plot([ang, ang], [1.05, 4.72], color="#edf3f8", linewidth=0.9, zorder=0)
        ax.text(
            ang,
            4.95,
            str(year),
            fontsize=9.4,
            color="#5b6b7c",
            ha="center",
            va="center",
            rotation=np.degrees(np.pi / 2 - ang),
            rotation_mode="anchor",
        )

    max_count = counts["count"].max() if not counts.empty else 1
    for _, row in counts.iterrows():
        year = int(row["Year_num"])
        state = row["Evidence_State"]
        count = int(row["count"])
        angle = angle_lookup[year]
        radius = radius_lookup[state]
        size = 40 + (count / max_count) * 2400
        ax.scatter(
            [angle],
            [radius],
            s=size,
            color=STATE_COLORS[state],
            alpha=0.90,
            edgecolors="white",
            linewidths=1.4,
            zorder=4,
        )
        if count >= max(4, int(np.ceil(max_count * 0.25))):
            ax.text(
                angle,
                radius,
                str(count),
                ha="center",
                va="center",
                fontsize=10.2,
                color="#16324f",
                fontweight="bold",
                zorder=5,
            )

    # Highlight the strongest evidence trajectory
    both_counts = (
        counts[counts["Evidence_State"] == "qPCR + RNA-seq"]
        .set_index("Year_num")["count"]
        .reindex(years, fill_value=0)
    )
    both_angles = [angle_lookup[y] for y in years if both_counts.loc[y] > 0]
    both_r = [radius_lookup["qPCR + RNA-seq"] for y in years if both_counts.loc[y] > 0]
    if both_angles:
        ax.plot(both_angles, both_r, color="#0f4c5c", linewidth=2.0, alpha=0.45, zorder=3)

    ax.text(
        0,
        0.22,
        "Validation\nEvidence\nChronomap",
        ha="center",
        va="center",
        fontsize=22,
        color="#16324f",
        fontweight="bold",
        linespacing=1.1,
    )
    ax.text(
        0,
        0.84,
        "Bubble size = studies",
        ha="center",
        va="center",
        fontsize=11.2,
        color="#6a7b8c",
    )

    grouped_labels = [
        ("Only in silico", 1.06, 0.70),
        ("qPCR only", 1.06, 0.62),
        ("RNA-seq only", 1.06, 0.54),
        ("qPCR + RNA-seq", 1.06, 0.46),
    ]
    for state, x, y in grouped_labels:
        ax.scatter(
            [x - 0.04],
            [y],
            s=170,
            color=STATE_COLORS[state],
            edgecolors="white",
            linewidths=1.2,
            transform=ax.transAxes,
            clip_on=False,
            zorder=6,
        )
        ax.text(
            x,
            y,
            state,
            transform=ax.transAxes,
            ha="left",
            va="center",
            fontsize=11.8,
            color="#16324f",
            fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.16", fc="white", ec="none", alpha=0.92),
            zorder=6,
        )

    fig.text(0.06, 0.955, "Validation and Expression Chronomap", fontsize=24, fontweight="bold", color="#16324f")

    plt.savefig(OUT_FILE, dpi=600, facecolor="#f4f7fb", bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)

    print("Validation chronomap generated.")
    print("Output:", OUT_FILE)


if __name__ == "__main__":
    main()
