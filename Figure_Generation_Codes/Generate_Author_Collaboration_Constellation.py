import os
import sys
import pandas as pd
import numpy as np
from collections import Counter
import itertools
import math
from repo_paths import configure_script_environment, resolve_data_file, resolve_output_dir

configure_script_environment()

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch, Circle
from matplotlib.path import Path
from matplotlib.colors import LinearSegmentedColormap

DATA_FILE = resolve_data_file("Data22.csv", "Data_cleaned.csv", "Data.csv")
OUT_DIR = resolve_output_dir()
OUT_FILE = os.path.join(OUT_DIR, "F8_Author_Collaboration_Constellation.png")

def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)

def main():
    ensure_dir(OUT_DIR)
    df = pd.read_csv(DATA_FILE, encoding='latin-1')

    authors_all = []
    pairs_all = []

    for idx, row in df.iterrows():
        if pd.isna(row['Authors']): continue
        raw_authors = str(row['Authors']).split(',')
        authors = [a.strip() for a in raw_authors if a.strip()]
        authors_all.extend(authors)
        if len(authors) > 1:
            for pair in itertools.combinations(sorted(authors), 2):
                pairs_all.append(pair)

    author_counts = Counter(authors_all)
    # Take top 35 authors for a beautifully dense, rich network
    top_authors = [item[0] for item in author_counts.most_common(35)]
    top_dict = {a: author_counts[a] for a in top_authors}

    pair_counts = Counter(pairs_all)
    # Filter pairs to only those where both authors are in our top 35 ring
    filtered_pairs = {p: c for p, c in pair_counts.items() if p[0] in top_authors and p[1] in top_authors}

    # Setup the Aesthetic Figure (White background)
    fig = plt.figure(figsize=(16, 16), facecolor="white")
    ax = fig.add_subplot(111, facecolor="white")
    ax.set_aspect('equal')
    ax.axis('off')

    # Calculate angles on a circle
    R = 10  # Radius
    n_nodes = len(top_authors)
    angles = np.linspace(0, 2 * np.pi, n_nodes, endpoint=False)
    
    node_coords = {}
    for i, author in enumerate(top_authors):
        x = R * math.cos(angles[i])
        y = R * math.sin(angles[i])
        node_coords[author] = (x, y, angles[i])

    # Draw Ribbons (Co-authorship connections)
    # Colormap for the ribbons based on intensity
    max_pair = max(filtered_pairs.values()) if filtered_pairs else 1
    
    for pair, count in filtered_pairs.items():
        a1, a2 = pair
        x1, y1, _ = node_coords[a1]
        x2, y2, _ = node_coords[a2]
        
        # Bezier curve control points
        # To make it organic, pull the curve toward the center (0,0)
        dist = math.hypot(x2 - x1, y2 - y1)
        # Closer nodes curve less, distant nodes curve deep into the center
        ctrl_pull = 0.5 * (1 - dist / (2*R)) if dist < (2*R) else 0
        cx, cy = ctrl_pull * (x1+x2)/2, ctrl_pull * (y1+y2)/2

        verts = [(x1, y1), (cx, cy), (x2, y2)]
        codes = [Path.MOVETO, Path.CURVE3, Path.CURVE3]
        
        # Colors: Deep contrasting blue adjusting by alpha
        intensity = min(1.0, count / max_pair)
        color = "#16324f"
        
        alpha = 0.15 + 0.6 * intensity
        lw = 1.0 + 4.0 * intensity

        patch = PathPatch(Path(verts, codes), facecolor='none', edgecolor=color, 
                          alpha=alpha, lw=lw, zorder=2, capstyle="round")
        ax.add_patch(patch)

    # Draw Nodes & Labels
    max_count = max(top_dict.values())
    for author, (x, y, ang) in node_coords.items():
        count = top_dict[author]
        intensity = count / max_count
        
        # Node Size
        radius = 0.2 + 0.6 * intensity
        
        # Inner Core
        ax.add_patch(Circle((x, y), radius, facecolor="#d1495b", edgecolor="white", lw=1.5, zorder=5))
        # Outer Ring
        ax.add_patch(Circle((x, y), radius*1.8, facecolor="#d1495b", alpha=0.15, zorder=4))
        
        # Text Labels Setup
        deg = math.degrees(ang)
        if 90 < deg < 270:
            rot = deg + 180
            ha = 'right'
            txt_x = x + (radius + 0.4) * math.cos(ang)
            txt_y = y + (radius + 0.4) * math.sin(ang)
        else:
            rot = deg
            ha = 'left'
            txt_x = x + (radius + 0.4) * math.cos(ang)
            txt_y = y + (radius + 0.4) * math.sin(ang)
        
        ax.text(txt_x, txt_y, f"{author} ({count})", 
                rotation=rot, rotation_mode="anchor",
                ha=ha, va="center", color="#16324f", fontsize=11, fontweight="bold", zorder=6)

    # Center Titles removed as requested

    # Bounding Box adjustments
    border = R + 4
    ax.set_xlim(-border, border)
    ax.set_ylim(-border, border)

    plt.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02)
    plt.savefig(OUT_FILE, dpi=400, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print(f"Author Constellation beautifully generated at: {OUT_FILE}")

if __name__ == "__main__":
    main()
