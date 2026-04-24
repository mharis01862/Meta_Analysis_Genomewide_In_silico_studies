import os
import re
import sys

from repo_paths import configure_script_environment, resolve_data_file, resolve_output_dir

configure_script_environment()
DATA_FILE = resolve_data_file("Data.csv", "Data_cleaned.csv", "Data22.csv")
OUT_RADIAL = os.path.join(resolve_output_dir(), "Innovative_Phylo_Astrolabe.png")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

COLOR_LOOKUP = {
    "Group":    "#f2c14e",
    "Subfamily": "#2f8393",
    "Clade":    "#49a7a3",
    "Class":    "#d1495b",
    "Subgroup": "#6c5ce7",
    "Type":     "#7d8597",
    "Cluster":  "#e17055",
    "Other":    "#97a6ba",
}

def clean_text(value):
    if pd.isna(value): return ""
    text = str(value).replace("\x96", "-").replace("\xa0", " ")
    return " ".join(text.strip().split())

def word_to_num(text):
    words = {
        "one":"1","two":"2","three":"3","four":"4","five":"5",
        "six":"6","seven":"7","eight":"8","nine":"9","ten":"10",
        "eleven":"11","twelve":"12","thirteen":"13","fourteen":"14",
        "fifteen":"15","sixteen":"16","seventeen":"17","eighteen":"18",
        "nineteen":"19","twenty":"20","twenty-one":"21","twenty-two":"22",
        "twenty-three":"23","twenty-four":"24","twenty-five":"25"
    }
    for w, n in words.items():
        text = re.sub(rf"\b{w}\b", n, text, flags=re.IGNORECASE)
    return text

def parse_gene_count(value):
    match = re.search(r"\d+(?:\.\d+)?", clean_text(value))
    return float(match.group(0)) if match else np.nan

def detect_structure_type(text):
    lower = clean_text(text).lower().replace("-","")
    if "subfamil" in lower: return "Subfamily"
    if "clade"    in lower: return "Clade"
    if "cluster"  in lower: return "Cluster"
    if "subgroup" in lower: return "Subgroup"
    if "group"    in lower: return "Group"
    if "class"    in lower: return "Class"
    if "type"     in lower: return "Type"
    if "not specified" in lower or lower == "": return "Not specified"
    return "Other"

def parse_partition_count(text):
    lower = word_to_num(clean_text(text).lower())
    if "not specified" in lower or lower == "": return np.nan
    patterns = [
        r"(?:into|in|classified|divided)\s*(?:into)?\s+(\d+)\s+(?:distinct\s+)?(?:main\s+|major\s+)?(?:subfamil|clade|group|class|subgroup|type|cluster|famil|lineage|clan)",
        r"^\s*(\d+)\s+(?:distinct\s+)?(?:main\s+|major\s+)?(?:subfamil|clade|group|class|subgroup|type|cluster|famil|lineage|clan)",
    ]
    for pat in patterns:
        m = re.search(pat, lower)
        if m: return float(m.group(1))
    for num in re.findall(r"\b(\d+)\b", lower):
        n = int(num)
        if 2 <= n <= 60: return float(n)
    return np.nan

def get_data():
    df = pd.read_csv(DATA_FILE, encoding="latin-1")
    df["Phylo_Text"]     = df["Phylogenetic Structure"].map(clean_text)
    df["Genes_num"]      = df["Number of Genes Identified"].map(parse_gene_count)
    df["Structure_Type"] = df["Phylo_Text"].map(detect_structure_type)
    df["Partition_Count"]= df["Phylo_Text"].map(parse_partition_count)
    plot_df = df[df["Structure_Type"] != "Not specified"].dropna(subset=["Partition_Count"]).copy()
    order = plot_df.groupby("Structure_Type").size().sort_values(ascending=False).index.tolist()
    plot_df["Structure_Type"] = pd.Categorical(plot_df["Structure_Type"], categories=order, ordered=True)
    return plot_df, order

def draw_radial_astrolabe(plot_df, order):
    fig = plt.figure(figsize=(15, 15), facecolor="white")
    ax  = fig.add_subplot(111, polar=True)
    ax.set_facecolor("white")

    n_cat = len(order)
    cat_angles  = np.linspace(0, 2*np.pi, n_cat, endpoint=False)
    angle_width = (2*np.pi) / n_cat
    max_p = float(plot_df["Partition_Count"].max())

    rng = np.random.default_rng(42)

    for i, typ in enumerate(order):
        color  = COLOR_LOOKUP.get(typ, "#97a6ba")
        subset = plot_df[plot_df["Structure_Type"] == typ]
        r_vals = subset["Partition_Count"].values

        ax.bar(cat_angles[i], max_p * 1.06, width=angle_width*0.82,
               color=color, alpha=0.15, bottom=0, zorder=0)

        jitter = rng.uniform(-angle_width*0.34, angle_width*0.34, size=len(r_vals))
        gene_v = subset["Genes_num"].fillna(subset["Genes_num"].median()).replace(0, 1)
        sizes  = np.clip(np.sqrt(gene_v) * 5.2, 30, 300)

        # Removed median circles. Increased alpha for visibility
        ax.scatter(cat_angles[i] + jitter, r_vals, s=sizes,
                   c=color, alpha=0.92, edgecolors=color, linewidths=1.0, zorder=3)

        label_r = max_p * 1.15
        
        # Merge text to avoid overlap between group name and count
        label_text = f"{typ.upper()}\n(n={len(subset)})"
        
        ax.text(cat_angles[i], label_r, label_text, ha="center", va="center",
                color=color, fontweight="bold", fontsize=14, zorder=4)

    tick_vals = [v for v in [5, 10, 15, 20, 30, 40, 50, 60] if v <= max_p]
    ax.set_yticks(tick_vals)
    ax.set_yticklabels([str(v) for v in tick_vals], color="#a0aabf", fontsize=10)
    ax.set_ylim(0, max_p * 1.25)
    ax.set_xticks([])
    ax.spines["polar"].set_visible(False)
    ax.grid(color="#dde6f2", linestyle="--", linewidth=1.2, zorder=1)

    fig.text(0.5, 0.95, "Phylogenetic Structure — Astrolabe Scatter",
             ha="center", va="top", fontsize=24, fontweight="bold", color="#16324f")
    fig.text(0.5, 0.91, "Each dot = one study  |  Dot size is proportional to genes identified",
             ha="center", va="top", fontsize=12, color="#6b7b8b")

    plt.savefig(OUT_RADIAL, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print("Astrolabe saved:", OUT_RADIAL)

def main():
    plot_df, order = get_data()
    print(f"Plotting {len(plot_df)} studies.")
    draw_radial_astrolabe(plot_df, order)

if __name__ == "__main__":
    main()
