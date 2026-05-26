"""
01_temporal_signal.py
---------------------
Temporal signal assessment via root-to-tip regression.

Reads the MCC nexus tree from code/, extracts collection year from tip
labels (format: ACCESSION_YYYY_MM_DD), computes root-to-tip distances,
runs linear regression against collection year, and saves:
  - results/temporal_signal_regression.csv   (per-tip data)
  - results/temporal_signal_plot.png         (scatter + regression line)
  - results/temporal_signal_summary.txt      (r, R², p-value, slope)

Usage:
    python code/01_temporal_signal.py --tree code/<your_tree>.tree
"""

import argparse
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo
from scipy import stats

RESULTS = Path("results")
RESULTS.mkdir(exist_ok=True)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--tree",
        required=True,
        help="Path to MCC nexus tree file, e.g. code/config2_fixed.tree",
    )
    return p.parse_args()


def extract_year(name: str) -> float | None:
    """Parse decimal year from tip label ACCESSION_YYYY_MM_DD or ACCESSION_YYYY."""
    m = re.search(r"_(\d{4})_(\d{2})_(\d{2})$", name)
    if m:
        y, mo, d = int(m.group(1)), int(m.group(2)), int(m.group(3))
        return y + (mo - 1) / 12 + (d - 1) / 365
    m = re.search(r"_(\d{4})$", name)
    if m:
        return float(m.group(1)) + 0.5  # mid-year imputation
    return None


def root_to_tip_distances(tree) -> dict[str, float]:
    """Return {tip_name: root-to-tip distance} for all terminals."""
    distances = {}
    for tip in tree.get_terminals():
        dist = tree.distance(tree.root, tip)
        distances[tip.name] = dist
    return distances


def main():
    args = parse_args()
    tree_path = Path(args.tree)
    if not tree_path.exists():
        sys.exit(f"Tree file not found: {tree_path}")

    print(f"Loading tree: {tree_path}")
    tree = Phylo.read(str(tree_path), "nexus")
    tree.root_at_midpoint()

    tips = tree.get_terminals()
    print(f"  {len(tips)} tips found")

    rtd = root_to_tip_distances(tree)

    rows = []
    skipped = 0
    for tip in tips:
        year = extract_year(tip.name)
        if year is None:
            skipped += 1
            continue
        rows.append({"tip": tip.name, "year": year, "root_to_tip": rtd[tip.name]})

    if skipped:
        print(f"  Skipped {skipped} tips with unparseable dates")

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "temporal_signal_regression.csv", index=False)

    # Linear regression
    slope, intercept, r, p, se = stats.linregress(df["year"], df["root_to_tip"])
    r2 = r**2

    summary = (
        f"Temporal signal summary\n"
        f"=======================\n"
        f"Tips used       : {len(df)}\n"
        f"Year range      : {df['year'].min():.1f} – {df['year'].max():.1f}\n"
        f"Slope           : {slope:.6e}  (subs/site/year)\n"
        f"Intercept       : {intercept:.6f}\n"
        f"Pearson r       : {r:.4f}\n"
        f"R²              : {r2:.4f}\n"
        f"p-value         : {p:.4e}\n"
        f"Std error       : {se:.6e}\n\n"
        f"Interpretation  :\n"
    )
    if r2 >= 0.3 and p < 0.05:
        summary += "  Good temporal signal — clock rate is likely data-driven.\n"
    elif p < 0.05:
        summary += "  Weak but significant signal — interpret TMRCA with caution.\n"
    else:
        summary += "  No significant temporal signal — BEAST estimates may be prior-dominated.\n"

    print(summary)
    (RESULTS / "temporal_signal_summary.txt").write_text(summary)

    # Plot
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(df["year"], df["root_to_tip"], s=12, alpha=0.5, color="#5B7EC9", label="Tips")
    x_line = np.linspace(df["year"].min(), df["year"].max(), 200)
    ax.plot(x_line, slope * x_line + intercept, color="#D85A30", linewidth=1.8,
            label=f"Regression (r={r:.3f}, R²={r2:.3f})")
    ax.set_xlabel("Collection year")
    ax.set_ylabel("Root-to-tip distance (subs/site)")
    ax.set_title("Temporal signal — HPV-16 MCC tree")
    ax.legend(fontsize=9)
    ax.text(0.03, 0.95, f"slope = {slope:.2e} subs/site/yr\np = {p:.3e}",
            transform=ax.transAxes, fontsize=8, va="top",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.7))
    fig.tight_layout()
    fig.savefig(RESULTS / "temporal_signal_plot.png", dpi=150)
    print("Saved results/temporal_signal_regression.csv")
    print("Saved results/temporal_signal_plot.png")
    print("Saved results/temporal_signal_summary.txt")


if __name__ == "__main__":
    main()
