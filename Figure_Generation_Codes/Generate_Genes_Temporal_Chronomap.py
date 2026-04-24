import os
import sys
import pandas as pd
import numpy as np
import re
from repo_paths import configure_script_environment, resolve_data_file, resolve_output_dir

configure_script_environment()

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import LinearSegmentedColormap

DATA_FILE = resolve_data_file("Data22.csv", "Data_cleaned.csv", "Data.csv")
OUT_DIR = resolve_output_dir()
OUT_FILE = os.path.join(OUT_DIR, "F10_Genes_Temporal_Chronomap.png")

def clean_text(value):
    return '' if pd.isna(value) else str(value).strip()

def parse_year(value):
    m = re.search(r'(19|20)\d{2}', clean_text(value))
    return int(m.group(0)) if m else None

def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)

def main():
    ensure_dir(OUT_DIR)
    df = pd.read_csv(DATA_FILE, encoding='latin-1')
    df['Year_num'] = df['Year'].map(parse_year)
    
    sub = df.dropna(subset=['Year_num', 'Number of Genes Identified']).copy()
    sub['Year_num'] = sub['Year_num'].astype(int)
    sub['Genes'] = pd.to_numeric(sub['Number of Genes Identified'], errors='coerce')
    sub = sub.dropna(subset=['Genes'])
    sub = sub[sub['Genes'] > 0] # ensure strictly positive for valid log scale
    
    # Sort logically
    sub = sub.sort_values("Year_num")

    fig, ax = plt.subplots(figsize=(15, 8.5), facecolor='white')
    ax.set_facecolor("white")
    
    # Jitter years horizontally for better bubble separation (prevents perfect overlap)
    np.random.seed(42)
    jitter = np.random.normal(0, 0.12, size=len(sub))
    x_data = sub['Year_num'] + jitter
    y_data = sub['Genes']
    
    # Premium deep-to-vibrant colormap exactly matching the established aesthetics
    aesthetic_cmap = LinearSegmentedColormap.from_list("aesthetic", ["#16324f", "#1f7a8c", "#49a5a1", "#ffd166", "#d1495b"])
    
    # Plot Bubbles (Size and Color scale dynamically with Log10 of the Gene Count)
    scatter = ax.scatter(x_data, y_data, 
                         s=15 + (np.log10(y_data) * 55), # Dynamic size scaling
                         c=np.log10(y_data), 
                         cmap=aesthetic_cmap,
                         alpha=0.75, 
                         edgecolors='white', 
                         linewidths=0.8,
                         zorder=4)

    # Calculate Trendline: Moving Average of the Median (resilient to severe outliers like 1900+)
    yearly_stats = sub.groupby('Year_num')['Genes'].median().reset_index()
    
    # Smooth the trendline using a rolling average
    yearly_stats['Genes_Smooth'] = yearly_stats['Genes'].rolling(window=2, min_periods=1, center=True).mean()
    
    # Plot the sweeping trendline
    ax.plot(yearly_stats['Year_num'], yearly_stats['Genes_Smooth'], color="#d1495b", linewidth=3.5, zorder=5, label="Moving Average Trajectory")
    
    # Fill area under the curve for visual anchoring
    ax.fill_between(yearly_stats['Year_num'], 0.5, yearly_stats['Genes_Smooth'], color="#d1495b", alpha=0.06, zorder=2)
    
    # Vertical Era Demarcation Lines matching the Results section text accurately
    ax.axvline(2012.5, color='#cbd5e1', linestyle='--', lw=2, zorder=3)
    ax.axvline(2016.5, color='#cbd5e1', linestyle='--', lw=2, zorder=3)
    ax.axvline(2020.5, color='#cbd5e1', linestyle='--', lw=2, zorder=3)
    
    # Era Labels positioned near the top of the chart
    y_label_pos = 1200
    ax.text(2010, y_label_pos, "Early Phase\n(Before 2013)", ha='center', color='#5b6b7c', fontweight='bold', fontsize=11, zorder=6)
    ax.text(2014.5, y_label_pos, "Expansion Phase\n(2013-2016)", ha='center', color='#5b6b7c', fontweight='bold', fontsize=11, zorder=6)
    ax.text(2018.5, y_label_pos, "Peak Growth\n(2017-2020)", ha='center', color='#5b6b7c', fontweight='bold', fontsize=11, zorder=6)
    ax.text(2023, y_label_pos, "Modern Era\n(2021-2025)", ha='center', color='#5b6b7c', fontweight='bold', fontsize=11, zorder=6)

    # Formatting Log Scale Y-Axis
    ax.set_yscale('log')
    ax.set_ylim(0.8, 3000)
    
    # Force Y-axis labels to standard readable integers instead of scientific 10^x
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f'{int(y)}'))
    
    # Clean up spines and grids
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('#cbd5e1')
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_color('#cbd5e1')
    ax.spines['bottom'].set_linewidth(1.5)
    ax.grid(axis='y', linestyle='-', alpha=0.3, color='#cbd5e1', zorder=1)

    # Title & Axis Labels
    ax.set_title("Temporal Chronomap of Genome-wide Gene Identifications", fontsize=20, fontweight='bold', color='#16324f', pad=25)
    ax.set_xlabel("Publication Year", fontsize=14, fontweight='500', color='#16324f', labelpad=15)
    ax.set_ylabel("Number of Genes Identified (Log Scale)", fontsize=14, fontweight='500', color='#16324f', labelpad=15)

    # Ensure X-axis only plots integer years cleanly
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax.tick_params(axis='x', colors='#16324f', labelsize=11)
    ax.tick_params(axis='y', colors='#16324f', labelsize=11)

    plt.tight_layout()
    plt.savefig(OUT_FILE, dpi=400, facecolor='white', bbox_inches='tight')
    plt.close()

    print(f"Genes Temporal Chronomap beautifully generated at: {OUT_FILE}")

if __name__ == '__main__':
    main()
