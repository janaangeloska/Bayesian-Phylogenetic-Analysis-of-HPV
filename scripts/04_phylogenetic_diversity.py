"""
04_phylogenetic_diversity.py
----------------------------
Computes phylogenetic diversity (Faith's PD) per country by extracting
country from tip labels or a metadata CSV, then summing branch lengths
in the minimum spanning subtree for each country's tips.

Also computes:
  - Mean Pairwise Distance (MPD) per country
  - Species richness (tip count) per country
  - PD per tip (Faith PD / n tips) — sampling-effort-corrected

Outputs:
  - results/phylo_diversity.csv         (per-country PD, MPD, n)
  - results/phylo_diversity_bar.png     (bar chart Faith PD by country)
  - results/phylo_diversity_scatter.png (PD vs tip count — sampling bias)

Usage:
    python code/04_phylogenetic_diversity.py --tree code/<your>.tree
    python code/04_phylogenetic_diversity.py --tree code/<your>.tree --meta code/metadata.csv
"""

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo

RESULTS = Path("results")
RESULTS.mkdir(exist_ok=True)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--tree", required=True)
    p.add_argument("--meta", default=None,
                   help="Optional CSV with 'accession' and 'country' columns")
    p.add_argument("--min_tips", type=int, default=3,
                   help="Minimum tips per country to include. Default: 3")
    return p.parse_args()


def accession_from_name(name: str) -> str:
    return name.split("_")[0]


def country_from_name(name: str) -> str | None:
    """Try to extract country from tip label if embedded after last pipe or underscore block."""
    # Header format: ACCESSION | organism | date | precision | country
    if "|" in name:
        parts = [p.strip() for p in name.split("|")]
        if len(parts) >= 5:
            return parts[4]
    return None


def load_country_map(meta_path: str) -> dict[str, str]:
    df = pd.read_csv(meta_path)
    df.columns = df.columns.str.lower()
    country_col = next((c for c in df.columns if "country" in c), None)
    acc_col = next((c for c in df.columns if "accession" in c), None)
    if not country_col or not acc_col:
        print("  Warning: could not find accession/country columns in metadata")
        return {}
    return dict(zip(df[acc_col].astype(str).str.strip(),
                    df[country_col].astype(str).str.strip()))


def subtree_branch_length(tree, tip_objects) -> float:
    """Faith's PD: sum of branch lengths in the minimal subtree spanning tips."""
    if len(tip_objects) == 1:
        return tree.distance(tree.root, tip_objects[0])
    mrca = tree.common_ancestor(tip_objects)
    total = 0.0
    for clade in mrca.find_clades():
        if clade.branch_length:
            total += clade.branch_length
    return total


def mean_pairwise_distance(tree, tip_objects) -> float:
    """Mean of all pairwise root-to-tip distances as a proxy for MPD."""
    if len(tip_objects) < 2:
        return 0.0
    dists = [tree.distance(tree.root, t) for t in tip_objects]
    n = len(dists)
    total = 0.0
    count = 0
    for i in range(n):
        for j in range(i + 1, n):
            total += abs(dists[i] - dists[j])
            count += 1
    return total / count if count else 0.0


def main():
    args = parse_args()
    tree_path = Path(args.tree)
    if not tree_path.exists():
        sys.exit(f"Tree not found: {tree_path}")

    print(f"Loading tree: {tree_path}")
    tree = Phylo.read(str(tree_path), "nexus")
    tips = tree.get_terminals()
    print(f"  {len(tips)} tips")

    country_map = {}
    if args.meta:
        print(f"Loading metadata: {args.meta}")
        country_map = load_country_map(args.meta)

    # Assign country per tip
    tip_country: dict[str, str] = {}
    for tip in tips:
        acc = accession_from_name(tip.name)
        if acc in country_map:
            tip_country[tip.name] = country_map[acc]
        else:
            c = country_from_name(tip.name)
            tip_country[tip.name] = c if c else "Unknown"

    # Group tips by country
    country_tips: dict[str, list] = defaultdict(list)
    for tip in tips:
        country_tips[tip_country[tip.name]].append(tip)

    print(f"  {len(country_tips)} countries/regions detected")

    rows = []
    for country, ctips in sorted(country_tips.items()):
        if len(ctips) < args.min_tips:
            continue
        pd_val = subtree_branch_length(tree, ctips)
        mpd = mean_pairwise_distance(tree, ctips)
        rows.append({
            "country": country,
            "n_tips": len(ctips),
            "faith_pd": pd_val,
            "pd_per_tip": pd_val / len(ctips),
            "mpd": mpd,
        })

    df = pd.DataFrame(rows).sort_values("faith_pd", ascending=False)
    df.to_csv(RESULTS / "phylo_diversity.csv", index=False)

    print(f"\nTop 10 countries by Faith PD:")
    print(df.head(10).to_string(index=False))

    # Bar chart: Faith PD by country (top 20)
    top = df.head(20)
    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.bar(range(len(top)), top["faith_pd"], color="#378ADD", alpha=0.8)
    ax2 = ax.twinx()
    ax2.plot(range(len(top)), top["n_tips"], "o--", color="#D85A30",
             markersize=5, linewidth=1, label="Tip count")
    ax.set_xticks(range(len(top)))
    ax.set_xticklabels(top["country"], rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Faith's PD (sum of branch lengths)")
    ax2.set_ylabel("Number of tips", color="#D85A30")
    ax2.tick_params(axis="y", labelcolor="#D85A30")
    ax.set_title("Phylogenetic diversity per country (Faith PD) — HPV-16")
    ax2.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(RESULTS / "phylo_diversity_bar.png", dpi=150)

    # Scatter: PD vs n_tips (reveals sampling bias)
    fig2, ax3 = plt.subplots(figsize=(7, 5))
    ax3.scatter(df["n_tips"], df["faith_pd"], s=40, alpha=0.7, color="#1D9E75")
    for _, row in df.head(15).iterrows():
        ax3.annotate(row["country"], (row["n_tips"], row["faith_pd"]),
                     fontsize=7, alpha=0.8, xytext=(3, 3), textcoords="offset points")
    ax3.set_xlabel("Number of tips (sampling effort)")
    ax3.set_ylabel("Faith's PD")
    ax3.set_title("PD vs sampling effort — sampling bias assessment")
    # Fit line
    if len(df) > 3:
        from scipy import stats
        slope, intercept, r, *_ = stats.linregress(df["n_tips"], df["faith_pd"])
        x = np.linspace(df["n_tips"].min(), df["n_tips"].max(), 100)
        ax3.plot(x, slope * x + intercept, "--", color="#888780",
                 linewidth=1, label=f"r={r:.2f}")
        ax3.legend(fontsize=8)
    fig2.tight_layout()
    fig2.savefig(RESULTS / "phylo_diversity_scatter.png", dpi=150)

    print("Saved results/phylo_diversity.csv")
    print("Saved results/phylo_diversity_bar.png")
    print("Saved results/phylo_diversity_scatter.png")


if __name__ == "__main__":
    main()
