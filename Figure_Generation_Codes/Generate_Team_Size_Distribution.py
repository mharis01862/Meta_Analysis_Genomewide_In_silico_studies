import os
import sys

from repo_paths import configure_script_environment, resolve_output_dir

configure_script_environment()

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT_DIR = resolve_output_dir()
OUT_FILE = os.path.join(OUT_DIR, "F9_Team_Size_Distribution.png")

def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)

def main():
    ensure_dir(OUT_DIR)
    
    # Data mapped precisely from the user's results section request
    # Categories ordered logically by team scope: Small to Macro
    categories = ['1 to 3 authors\n(Small teams)', 
                  '4 to 6 authors\n(Moderate teams)', 
                  '7 to 10 authors\n(Standard assembly)', 
                  '> 10 authors\n(Macro-teams)']
    counts = [35, 106, 126, 42]
    
    # Very light pastel academic color palette
    colors = ["#cbd5e8", "#b3e2cd", "#fdcdac", "#f4cae4"]

    fig, ax = plt.subplots(figsize=(10.5, 8.5), facecolor="white")
    ax.set_facecolor("white")

    wedges, texts, autotexts = ax.pie(
        counts, 
        labels=categories, 
        colors=colors,
        autopct=lambda pct: f"{pct:.1f}%\n(n={int(round(pct/100.*sum(counts)))})",
        startangle=160,
        pctdistance=0.74,
        labeldistance=1.1,
        wedgeprops=dict(width=0.48, edgecolor='white', linewidth=4),
        textprops=dict(color="#16324f", fontsize=12, fontweight="bold")
    )
    
    # Update text inside the donut conditionally for contrast reading
    for i, autotext in enumerate(autotexts):
        autotext.set_color('#16324f')
        autotext.set_fontsize(13)
        autotext.set_fontweight('bold')

    # Titles removed per user request

    # Adding a sleek inner-circle annotation highlighting the total and median
    ax.text(0, 0, "Average Team\n7.03 Authors", ha='center', va='center', fontsize=18, fontweight='bold', color="#118ab2")

    plt.tight_layout(pad=3.0)
    plt.savefig(OUT_FILE, dpi=400, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print(f"Team Size Distribution cleanly generated at: {OUT_FILE}")

if __name__ == "__main__":
    main()
