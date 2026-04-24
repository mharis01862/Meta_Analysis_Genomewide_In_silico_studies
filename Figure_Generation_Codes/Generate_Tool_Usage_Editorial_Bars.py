import os
import re
import sys
from collections import Counter

from repo_paths import configure_script_environment, resolve_data_file, resolve_output_dir

configure_script_environment()
DATA_FILE = resolve_data_file("Data.csv", "Data_cleaned.csv", "Data22.csv")
OUT_FILE = os.path.join(resolve_output_dir(), "F8_tool_usage_editorial_bars.png")

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd


TOP_N_METHODS = 8

GROUP_COLORS = {
    "Sequence Identification": "#1f7a8c",
    "Domain Analysis": "#2a9d8f",
    "Motif Analysis": "#f4a261",
    "Phylogeny": "#3d5a98",
    "Synteny": "#8d5fd3",
    "Visualization": "#d1495b",
    "Other Methodology Tools": "#5f6b76",
}

GROUP_BG = {
    "Sequence Identification": "#d8e7ea",
    "Domain Analysis": "#d9eee9",
    "Motif Analysis": "#f8e6d4",
    "Phylogeny": "#dde4f1",
    "Synteny": "#e7def5",
    "Visualization": "#f3dce1",
    "Other Methodology Tools": "#e3e6e8",
}


def clean(x):
    if pd.isna(x):
        return ""
    return " ".join(str(x).strip().split())


def yes(x):
    t = clean(x).lower()
    return t not in {"", "no", "none", "false", "n/a", "na"} and "not used" not in t


def tokens(x):
    x = clean(x)
    if not x or x.lower() in {"no", "none", "not specified", "not available"}:
        return []
    x = re.sub(r"\([^)]*\)", "", x)
    x = x.replace(" and ", ", ")
    out = []
    for part in re.split(r"[;,/]", x):
        part = clean(part).strip(".")
        if len(part) >= 2 and part.lower() not in {"yes", "no", "none", "na", "n/a"}:
            out.append(part)
    return list(dict.fromkeys(out))


def top_tool(df, col):
    counter = Counter()
    for _, row in df.iterrows():
        for tool in set(tokens(row.get(col, ""))):
            counter[tool] += 1
    return counter.most_common(1)[0] if counter else ("None", 0)


def collect_rows():
    df = pd.read_csv(DATA_FILE, encoding="latin-1")

    hmm_vals = df["HMMER (with Pfam domains)"].fillna("").map(clean)
    blast_vals = df["BLAST (E-value cut-off used)"].fillna("").map(clean)

    hmm_total = sum(hmm_vals.map(yes))
    hmm_pfam_spec = 143  # Manually validated count of studies specifically providing the PFAM ID
    blast_total = sum(blast_vals.map(yes))
    blast_eval_spec = sum(
        1
        for v in blast_vals
        if yes(v)
        and (
            re.search(r"e-?value", v, re.I)
            or re.search(r"\b1e-\d+|\d+e-\d+|10-\d+", v, re.I)
            or re.search(r"0\.0*\d", v)
        )
    )

    other_domain_name, other_domain_count = top_tool(df, "Other (Tools for Domain Analysis)")
    other_motif_name, other_motif_count = top_tool(df, "Others (Tools for Motif Analysis)")
    other_synteny_name, other_synteny_count = top_tool(df, "Other (Tools for Synteny)")
    other_visual_name, other_visual_count = top_tool(df, "Other (Tools for Visualization)")

    rows = [
        ("Sequence Identification", "HMMER", hmm_total),
        ("Sequence Identification", "HMMER with Pfam/domain specified", hmm_pfam_spec),
        ("Sequence Identification", "BLAST", blast_total),
        ("Sequence Identification", "BLAST with E-value specified", blast_eval_spec),
        ("Domain Analysis", "SMART", sum(df["SMART"].map(yes))),
        ("Domain Analysis", "CDD", sum(df["CDD"].map(yes))),
        ("Domain Analysis", f"Other domain tool: {other_domain_name}", other_domain_count),
        ("Motif Analysis", "MEME", sum(df["MEME"].map(yes))),
        ("Motif Analysis", f"Other motif tool: {other_motif_name}", other_motif_count),
        ("Phylogeny", "MEGA", sum(df["MEGA"].map(yes))),
        ("Phylogeny", "IQ-TREE", sum(df["IQ-TREE"].map(yes))),
        ("Synteny", "MCScanX", sum(df["MCScanX"].map(yes))),
        ("Synteny", f"Other synteny tool: {other_synteny_name}", other_synteny_count),
        ("Visualization", "TBTools", sum(df["TBTools"].map(yes))),
        ("Visualization", f"Other visualization tool: {other_visual_name}", other_visual_count),
    ]

    methodology_counter = Counter()
    methodology_col = "Other Tools Used In Methodology"
    already_shown = {
        "mega", "tbtools", "pfam", "gsds", "meme", "mcscanx",
        "smart", "cdd", "blast", "hmmer", "iq-tree", "weblogo",
        "plant genome duplication database",
    }
    db_keywords = {"database", "plantcare", "pfam", "string", "interpro", "interproscan", "pgdd"}
    for _, row in df.iterrows():
        for tool in set(tokens(row.get(methodology_col, ""))):
            low = tool.lower()
            if any(k == low or k in low for k in db_keywords):
                continue
            if low in already_shown:
                continue
            methodology_counter[tool] += 1
    for tool, count in methodology_counter.most_common(TOP_N_METHODS):
        rows.append(("Other Methodology Tools", tool, count))

    return rows


def main():
    rows = collect_rows()
    group_order = [
        "Sequence Identification",
        "Domain Analysis",
        "Motif Analysis",
        "Phylogeny",
        "Synteny",
        "Visualization",
        "Other Methodology Tools",
    ]

    ordered_rows = []
    for group in group_order:
        group_rows = [r for r in rows if r[0] == group]
        group_rows.sort(key=lambda r: (-r[2], r[1].lower()))
        if group_rows:
            ordered_rows.extend(group_rows)

    n = len(ordered_rows)
    max_count = max(r[2] for r in ordered_rows)
    positions = np.arange(n)[::-1]

    fig_h = max(10.5, 0.5 * n + 2.0)
    fig, ax = plt.subplots(figsize=(16.4, fig_h), facecolor="white")
    ax.set_facecolor("white")

    group_bounds = {}
    for group in group_order:
        idxs = [i for i, r in enumerate(ordered_rows) if r[0] == group]
        if idxs:
            ys = positions[idxs]
            group_bounds[group] = (ys.max() + 0.6, ys.min() - 0.6)

    for group, (y_top, y_bottom) in group_bounds.items():
        color = GROUP_COLORS[group]
        bg = GROUP_BG[group]
        rect = patches.FancyBboxPatch(
            (-max_count * 0.03, y_bottom),
            max_count * 1.08,
            y_top - y_bottom,
            boxstyle="round,pad=0.02,rounding_size=0.22",
            linewidth=0,
            facecolor=bg,
            alpha=0.75,
            zorder=0,
        )
        ax.add_patch(rect)
        ax.text(
            -max_count * 0.17,
            (y_top + y_bottom) / 2,
            group,
            ha="left",
            va="center",
            fontsize=12.8,
            fontweight="bold",
            color="#16324f",
            bbox=dict(boxstyle="round,pad=0.20", fc="white", ec="none", alpha=0.96),
        )

    labels = []
    for (group, label, count), y in zip(ordered_rows, positions):
        color = GROUP_COLORS[group]
        ax.barh(y, count, color=color, height=0.58, alpha=0.92, zorder=2)
        ax.barh(y, count, color="none", edgecolor="white", linewidth=1.2, height=0.58, zorder=3)
        ax.scatter(count, y, s=90 + count * 4.2, color=color, edgecolors="white", linewidths=1.5, zorder=4)
        ax.text(
            count + max_count * 0.018,
            y,
            str(int(count)),
            ha="left",
            va="center",
            fontsize=11.0,
            color="#34495e",
            fontweight="bold",
            zorder=5,
        )
        labels.append(label)

    ax.set_yticks(positions)
    ax.set_yticklabels(labels, fontsize=11.4, color="#16324f")
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(-max_count * 0.22, max_count * 1.08)
    ax.set_xlabel("Number of studies using the tool", fontsize=13.8, color="#16324f")
    ax.grid(axis="x", color="#e3ebf3", linewidth=0.95, alpha=0.9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color("#d7e2ee")

    plt.tight_layout()
    plt.savefig(OUT_FILE, dpi=340, facecolor="white", bbox_inches="tight")
    plt.close(fig)

    print("Editorial bar tool figure generated.")
    print("Output:", OUT_FILE)


if __name__ == "__main__":
    main()
