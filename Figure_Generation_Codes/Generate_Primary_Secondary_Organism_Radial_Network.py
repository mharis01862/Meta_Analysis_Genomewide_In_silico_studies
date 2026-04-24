import csv
import math
import os
import re
from collections import Counter

import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".pydeps"))
os.environ["MPLCONFIGDIR"] = os.path.join(os.path.dirname(__file__), ".mplconfig")

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, PathPatch
from matplotlib.path import Path

from repo_paths import configure_script_environment, resolve_data_file, resolve_output_dir

configure_script_environment()
DATA_FILE = resolve_data_file("Data_cleaned.csv", "Data22.csv", "Data.csv")
OUT_DIR = resolve_output_dir()
OUT_FILE = os.path.join(OUT_DIR, "Primary_Secondary_Organism_Radial_Network.png")


PRIMARY_COLOR = "#0b4f8c"
PRIMARY_GLOW = "#8ecae6"
SECONDARY_COLOR = "#2a9d8f"
SECONDARY_GLOW = "#b7e4d6"
EDGE_COLOR = "#79aeda"
TEXT_DARK = "#16324f"
TEXT_MUTED = "#557089"


def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)


def clean_text(text):
    text = "" if text is None else str(text).strip()
    text = re.sub(r"\s+", " ", text)
    return text.strip(" ,;")


def compact(text, n=40):
    text = clean_text(text)
    return text if len(text) <= n else text[: n - 3] + "..."


def split_secondary(text):
    text = clean_text(text)
    if not text or text.lower() in {"none", "nan", "not specified"}:
        return []
    text = text.replace(" and ", ", ")
    parts = []
    for token in re.split(r"[;,/]", text):
        token = clean_text(token)
        if token and token.lower() not in {"none", "nan"}:
            parts.append(token)
    return parts


def load_network_data():
    with open(DATA_FILE, encoding="latin-1", newline="") as f:
        rows = list(csv.DictReader(f))

    primary_counts = Counter()
    secondary_counts = Counter()
    edge_counts = Counter()

    for row in rows:
        primary = clean_text(row.get("Primary Organism", ""))
        secondaries = split_secondary(row.get("Secondary Organism", ""))
        if primary:
            primary_counts[primary] += 1
        for sec in secondaries:
            secondary_counts[sec] += 1
            if primary:
                edge_counts[(primary, sec)] += 1

    return primary_counts, secondary_counts, edge_counts


def polar_positions(labels, radius, start_deg, end_deg):
    if len(labels) == 1:
        angles = [(start_deg + end_deg) / 2]
    else:
        step = (end_deg - start_deg) / (len(labels) - 1)
        angles = [start_deg + i * step for i in range(len(labels))]
    pos = {}
    for label, angle in zip(labels, angles):
        rad = math.radians(angle)
        pos[label] = (0.5 + radius * math.cos(rad), 0.5 + radius * math.sin(rad), angle)
    return pos


def draw_curve(ax, start, end, weight):
    x0, y0, _ = start
    x1, y1, _ = end
    ctrl1 = (x0 * 0.72 + 0.5 * 0.28, y0 * 0.72 + 0.5 * 0.28)
    ctrl2 = (x1 * 0.72 + 0.5 * 0.28, y1 * 0.72 + 0.5 * 0.28)
    path = Path(
        [(x0, y0), ctrl1, ctrl2, (x1, y1)],
        [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4],
    )
    glow = PathPatch(
        path,
        facecolor="none",
        edgecolor=PRIMARY_GLOW,
        linewidth=1.8 + weight * 0.9,
        alpha=0.04 + weight * 0.013,
        capstyle="round",
        zorder=1,
    )
    core = PathPatch(
        path,
        facecolor="none",
        edgecolor=EDGE_COLOR,
        linewidth=0.5 + weight * 0.34,
        alpha=0.10 + weight * 0.024,
        capstyle="round",
        zorder=2,
    )
    ax.add_patch(glow)
    ax.add_patch(core)


def draw_node(ax, x, y, radius, fill, glow, label, angle, count_value, role):
    ax.add_patch(Circle((x, y), radius=radius * 1.45, facecolor=glow, edgecolor="none", alpha=0.22, zorder=3))
    ax.add_patch(Circle((x, y), radius=radius, facecolor=fill, edgecolor="white", linewidth=1.4, zorder=4))

    if role == "primary":
        offset = 0.04
        ha = "left" if math.cos(math.radians(angle)) >= 0 else "right"
        tx = x + offset if ha == "left" else x - offset
        ty = y
    else:
        ha = "center"
        tx = x
        ty = y - radius - 0.032

    label_text = ax.text(tx, ty + 0.008, label, ha=ha, va="center", fontsize=8.5, color=TEXT_DARK, fontstyle="italic", fontweight="bold", zorder=5)
    count_text = ax.text(tx, ty - 0.016, f"n={count_value}", ha=ha, va="center", fontsize=8.0, color=TEXT_MUTED, zorder=5)
    return label_text, count_text


def main():
    ensure_dir(OUT_DIR)
    primary_counts, secondary_counts, edge_counts = load_network_data()

    top_primary = [k for k, _ in primary_counts.most_common(16)]
    top_secondary = [k for k, _ in secondary_counts.most_common(10)]
    filtered_edges = {
        (p, s): w
        for (p, s), w in edge_counts.items()
        if p in top_primary and s in top_secondary and w >= 3
    }

    used_primary = [p for p in top_primary if any((p, s) in filtered_edges for s in top_secondary)]
    used_secondary = [s for s in top_secondary if any((p, s) in filtered_edges for p in top_primary)]

    fig, ax = plt.subplots(figsize=(8, 7), facecolor="white")
    ax.set_facecolor("white")
    ax.set_aspect("equal")

    # Subtle guide rings
    ax.add_patch(Circle((0.5, 0.5), radius=0.30, facecolor="white", edgecolor="#e4eef7", linewidth=1.0, zorder=0))
    ax.add_patch(Circle((0.5, 0.5), radius=0.42, facecolor="none", edgecolor="#edf3f8", linewidth=1.0, zorder=0))

    primary_pos = polar_positions(used_primary, radius=0.42, start_deg=118, end_deg=422)
    secondary_pos = polar_positions(used_secondary, radius=0.19, start_deg=130, end_deg=410)

    secondary_label_overrides = {
        "Solanum lycopersicum": (0.0, -0.028),
        "Populus trichocarpa": (0.0, 0.03),
    }

    max_edge = max(filtered_edges.values()) if filtered_edges else 1
    for (p, s), weight in sorted(filtered_edges.items(), key=lambda kv: kv[1]):
        draw_curve(ax, primary_pos[p], secondary_pos[s], weight / max_edge * 8)

    max_primary = max(primary_counts[p] for p in used_primary) if used_primary else 1
    max_secondary = max(secondary_counts[s] for s in used_secondary) if used_secondary else 1

    for p in used_primary:
        x, y, angle = primary_pos[p]
        radius = 0.011 + 0.017 * (primary_counts[p] / max_primary)
        draw_node(ax, x, y, radius, PRIMARY_COLOR, PRIMARY_GLOW, compact(p, 24), angle, primary_counts[p], "primary")

    for s in used_secondary:
        x, y, angle = secondary_pos[s]
        radius = 0.012 + 0.022 * (secondary_counts[s] / max_secondary)
        label_text, count_text = draw_node(ax, x, y, radius, SECONDARY_COLOR, SECONDARY_GLOW, compact(s, 22), angle, secondary_counts[s], "secondary")
        if s in secondary_label_overrides:
            dx, dy = secondary_label_overrides[s]
            label_text.set_position((label_text.get_position()[0] + dx, label_text.get_position()[1] + dy))
            count_text.set_position((count_text.get_position()[0] + dx, count_text.get_position()[1] + dy))

    ax.text(0.5, 0.515, "Comparative\nReference\nCore", ha="center", va="center", fontsize=15, color=SECONDARY_COLOR, fontweight="bold")
    ax.text(
        0.5,
        0.735,
        "Secondary Reference Organisms",
        ha="center",
        va="center",
        fontsize=14,
        color=SECONDARY_COLOR,
        fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#d7efe7", alpha=0.95),
        zorder=6,
    )
    ax.text(
        0.5,
        0.045,
        "Primary Organisms",
        ha="center",
        va="center",
        fontsize=14,
        color=PRIMARY_COLOR,
        fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#dceaf6", alpha=0.95),
        zorder=6,
    )



    ax.set_title(
        "Primary-Secondary Organism Radial Comparison Network",
        loc="left",
        fontsize=22,
        pad=18,
        color=TEXT_DARK,
        fontweight="bold",
    )
    ax.text(
        0.0,
        -0.035,
        "Outer nodes are study organisms, inner nodes are recurrent reference organisms, and curved links represent repeated comparative pairing across studies.",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        color=TEXT_MUTED,
    )

    ax.set_xlim(-0.35, 1.35)
    ax.set_ylim(-0.30, 1.30)
    ax.axis("off")
    # Removed subplots_adjust as bbox_inches='tight' is used
    plt.savefig(OUT_FILE, dpi=150, facecolor="white", bbox_inches="tight", pad_inches=0.2)
    plt.close(fig)

    print("Primary-secondary radial network generated successfully.")
    print(f"OUTPUT DIRECTORY: {OUT_DIR}")
    print(f"OUTPUT FILE: {OUT_FILE}")


if __name__ == "__main__":
    main()
 








