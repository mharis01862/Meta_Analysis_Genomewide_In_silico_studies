import csv
import json
import os
import sys
from collections import Counter

from repo_paths import configure_script_environment, resolve_auxiliary_file, resolve_data_file, resolve_output_dir

configure_script_environment()

DATA_FILE = resolve_data_file("Data22.csv", "Data_cleaned.csv", "Data.csv")
GEOJSON_FILE = resolve_auxiliary_file("world_countries.geojson")
OUT_FILE = os.path.join(resolve_output_dir(), "F1_CountryOrigin.png")

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.colors import LinearSegmentedColormap, PowerNorm
from matplotlib.patches import Polygon


COUNTRY_MAP = {
    "United States": "United States of America",
    "USA": "United States of America",
    "US": "United States of America",
    "Korea": "South Korea",
}


def normalize_country(value):
    text = (value or "").strip()
    if not text or text.lower() == "unknown":
        return ""
    return COUNTRY_MAP.get(text, text)


def load_counts():
    with open(DATA_FILE, encoding="latin-1", newline="") as f:
        rows = list(csv.DictReader(f))
    counts = Counter()
    for row in rows:
        country = normalize_country(row.get("Country of Origin", ""))
        if country:
            counts[country] += 1
    return counts


def iter_polygons(geometry):
    gtype = geometry["type"]
    coords = geometry["coordinates"]
    if gtype == "Polygon":
        for ring in coords[:1]:
            yield ring
    elif gtype == "MultiPolygon":
        for poly in coords:
            for ring in poly[:1]:
                yield ring


def main():
    counts = load_counts()
    with open(GEOJSON_FILE, encoding="utf-8") as f:
        geo = json.load(f)

    cmap = LinearSegmentedColormap.from_list(
        "country_pub_visible",
        ["#eef8d8", "#cce98e", "#98d356", "#58ac2a", "#166728"],
    )
    norm = PowerNorm(gamma=0.42, vmin=1, vmax=max(counts.values()) if counts else 1)

    fig, ax = plt.subplots(figsize=(22, 10.8), facecolor="white")
    ax.set_facecolor("white")

    patches = []
    colors = []
    label_points = []

    for feature in geo["features"]:
        name = feature["properties"].get("name", "")
        count = counts.get(name, 0)
        face = "#f1f1f1" if count == 0 else cmap(norm(count))
        rings = list(iter_polygons(feature["geometry"]))
        for ring in rings:
            patches.append(Polygon(ring, closed=True))
            colors.append(face)
        if count > 0 and rings:
            largest = max(rings, key=lambda pts: len(pts))
            xs = [p[0] for p in largest]
            ys = [p[1] for p in largest]
            label_points.append((name, count, sum(xs) / len(xs), sum(ys) / len(ys)))

    pc = PatchCollection(
        patches,
        facecolor=colors,
        edgecolor="#d9d9d9",
        linewidths=0.6,
        zorder=2,
    )
    ax.add_collection(pc)

    nudges = {
        "China": (9, 2),
        "India": (8, -2),
        "Pakistan": (3, 3),
        "South Korea": (11, 2),
        "Bangladesh": (8, -1),
        "Turkey": (7, 3),
        "Brazil": (-10, -6),
        "Iran": (7, 0),
        "United States of America": (-19, 8),
        "Canada": (-21, 11),
        "Australia": (11, -2),
        "Mexico": (-10, -2),
        "Saudi Arabia": (10, -1),
        "Denmark": (6, 6),
        "Italy": (8, -2),
        "Austria": (8, 4),
        "Netherlands": (8, 6),
        "Spain": (-10, -2),
        "Belarus": (10, 6),
        "Malaysia": (9, -1),
        "United Arab Emirates": (10, 2),
        "Sudan": (8, -2),
        "Greece": (8, -2),
    }

    for name, count, x, y in sorted(label_points, key=lambda item: item[1], reverse=True):
        dx, dy = nudges.get(name, (4, 4))
        tx, ty = x + dx, y + dy
        ax.plot([x, tx], [y, ty], color="#7a8694", linewidth=0.65, zorder=3, alpha=0.8)
        ax.text(
            tx,
            ty,
            f"{name.replace('United States of America', 'USA')} ({count})",
            fontsize=10.1,
            color="#243746",
            ha="left" if dx >= 0 else "right",
            va="center",
            bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="#d7e6c5", alpha=0.96),
            zorder=4,
        )

    ax.set_xlim(-180, 180)
    ax.set_ylim(-60, 90)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.022, pad=0.018, shrink=0.78)
    cbar.set_label("Number of publications", fontsize=12.5)
    cbar.ax.tick_params(labelsize=11)

    plt.subplots_adjust(left=0.01, right=0.94, top=0.98, bottom=0.03)
    plt.savefig(OUT_FILE, dpi=350, bbox_inches="tight", facecolor="white", pad_inches=0.02)
    plt.close(fig)

    print("Country origin final map regenerated.")
    print("Output:", OUT_FILE)


if __name__ == "__main__":
    main()
