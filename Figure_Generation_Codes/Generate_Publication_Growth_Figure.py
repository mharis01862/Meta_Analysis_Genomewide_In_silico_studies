import os
import re
import sys

from repo_paths import configure_script_environment, resolve_data_file, resolve_output_dir

configure_script_environment()

DATA_FILE = resolve_data_file("Data.csv", "Data_cleaned.csv", "Data22.csv")
OUT_FILE = os.path.join(resolve_output_dir(), "F2_Yearwise_growth.png")

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd


plt.rcParams["figure.facecolor"] = "#ffffff"
plt.rcParams["axes.facecolor"] = "#ffffff"
plt.rcParams["savefig.facecolor"] = "#ffffff"
plt.rcParams["axes.edgecolor"] = "#d7e2ee"
plt.rcParams["grid.color"] = "#dfe7f1"
plt.rcParams["axes.titleweight"] = "bold"


def parse_year(value):
    match = re.search(r"(19|20)\d{2}", str(value))
    return int(match.group(0)) if match else None


def main():
    df = pd.read_csv(DATA_FILE, encoding="latin-1")
    df["Year_num"] = df["Year"].map(parse_year)

    yearly = (
        df.dropna(subset=["Year_num"])
        .groupby("Year_num")
        .size()
        .reset_index(name="papers")
        .sort_values("Year_num")
    )
    yearly["cumulative"] = yearly["papers"].cumsum()
    peak = yearly.loc[yearly["papers"].idxmax()]

    fig, ax1 = plt.subplots(figsize=(15.2, 7.6), facecolor="#ffffff")
    ax2 = ax1.twinx()
    ax1.set_facecolor("#ffffff")
    ax2.set_facecolor("#ffffff")

    ax1.bar(
        yearly["Year_num"],
        yearly["papers"],
        width=0.75,
        color="#9ecae1",
        edgecolor="#0b4f8c",
        linewidth=1.0,
        alpha=0.9,
    )
    ax1.plot(
        yearly["Year_num"],
        yearly["papers"],
        color="#0b4f8c",
        linewidth=2.8,
        marker="o",
        markersize=5,
    )
    ax2.plot(
        yearly["Year_num"],
        yearly["cumulative"],
        color="#d62839",
        linewidth=3,
        linestyle="--",
    )

    ax1.scatter(
        [peak["Year_num"]],
        [peak["papers"]],
        s=120,
        color="#d62839",
        edgecolor="white",
        linewidth=1.5,
        zorder=5,
    )
    ax1.annotate(
        f"Peak output: {int(peak['Year_num'])}\n{int(peak['papers'])} studies",
        xy=(peak["Year_num"], peak["papers"]),
        xytext=(-118, 14),
        textcoords="offset points",
        fontsize=11,
        color="#7f1d1d",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#f0c2c2", alpha=0.95),
    )

    ax1.set_xlabel("Publication year")
    ax1.set_ylabel("Studies published per year")
    ax2.set_ylabel("Cumulative studies")
    ax1.set_xticks(yearly["Year_num"].astype(int).tolist())
    ax1.set_xticklabels([str(int(x)) for x in yearly["Year_num"]], rotation=45, ha="right")

    ax1.spines["top"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    ax1.grid(axis="y", alpha=0.6)
    ax1.grid(axis="x", visible=False)

    plt.tight_layout()
    plt.savefig(OUT_FILE, dpi=300, bbox_inches="tight", facecolor="#ffffff")
    plt.close(fig)

    print("Publication growth figure regenerated.")
    print("Output:", OUT_FILE)


if __name__ == "__main__":
    main()
