"""
03_tmrca_by_clade.py
--------------------
Extracts TMRCA (node height) and 95% HPD intervals from BEAST-annotated
MCC nexus tree internal node comments, then reports per-clade divergence
times relative to the most recent sample year.

BEAST annotates internal nodes with comments like:
  [&height=315.5,height_95%_HPD={235.6,428.1},posterior=0.98,...]

Outputs:
  - results/tmrca_nodes.csv          (all annotated internal nodes)
  - results/tmrca_by_clade.csv       (one row per major clade)
  - results/tmrca_plot.png           (HPD interval bar chart)
  - results/tmrca_summary.txt        (human-readable summary)

Usage:
    python code/03_tmrca_by_clade.py --tree code/<your>.tree
    python code/03_tmrca_by_clade.py --tree code/<your>.tree --ref_year 2024
"""

import argparse
import re
import sys
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
    p.add_argument("--ref_year", type=float, default=2024.5,
                   help="Most recent sample year (decimal). Default: 2024.5")
    return p.parse_args()


def parse_comment(comment: str) -> dict:
    """Extract key=value pairs from a BEAST nexus node comment string."""
    if not comment:
        return {}
    data = {}
    # height
    m = re.search(r"height=([\d.eE+\-]+)", comment)
    if m:
        data["height"] = float(m.group(1))
    # height_95%_HPD={lo,hi}
    m = re.search(r"height_95%_HPD=\{([\d.eE+\-]+),([\d.eE+\-]+)\}", comment)
    if m:
        data["height_hpd_lo"] = float(m.group(1))
        data["height_hpd_hi"] = float(m.group(2))
    # posterior
    m = re.search(r"posterior=([\d.eE+\-]+)", comment)
    if m:
        data["posterior"] = float(m.group(1))
    return data


def height_to_calendar(height: float, ref_year: float) -> float:
    return ref_year - height


def main():
    args = parse_args()
    tree_path = Path(args.tree)
    if not tree_path.exists():
        sys.exit(f"Tree not found: {tree_path}")

    ref_year = args.ref_year
    print(f"Loading tree: {tree_path}")
    print(f"Reference year: {ref_year}")

    tree = Phylo.read(str(tree_path), "nexus")
    internals = tree.get_nonterminals()
    print(f"  {len(tree.get_terminals())} tips, {len(internals)} internal nodes")

    rows = []
    for i, node in enumerate(internals):
        comment = getattr(node, "comment", None) or ""
        parsed = parse_comment(comment)
        if not parsed:
            continue
        h = parsed.get("height")
        lo = parsed.get("height_hpd_lo")
        hi = parsed.get("height_hpd_hi")
        post = parsed.get("posterior")
        ntips = len(node.get_terminals())
        rows.append({
            "node_id": i,
            "n_tips": ntips,
            "height_yrs": h,
            "hpd_lo_yrs": lo,
            "hpd_hi_yrs": hi,
            "calendar_mean": height_to_calendar(h, ref_year) if h else None,
            "calendar_hpd_lo": height_to_calendar(hi, ref_year) if hi else None,
            "calendar_hpd_hi": height_to_calendar(lo, ref_year) if lo else None,
            "posterior": post,
        })

    df = pd.DataFrame(rows)
    df_sorted = df.sort_values("n_tips", ascending=False)
    df_sorted.to_csv(RESULTS / "tmrca_nodes.csv", index=False)
    print(f"  {len(df)} annotated internal nodes extracted")

    # Identify key clades by tip count thresholds
    total = len(tree.get_terminals())
    thresholds = {
        "Root (whole tree)": (total * 0.90, total),
        "Major clade (>50%)": (total * 0.50, total * 0.90),
        "Mid clade (20-50%)": (total * 0.20, total * 0.50),
        "Small clade (5-20%)": (total * 0.05, total * 0.20),
        "Mini clade (<5%)": (0, total * 0.05),
    }

    clade_rows = []
    for label, (lo_n, hi_n) in thresholds.items():
        subset = df[(df["n_tips"] > lo_n) & (df["n_tips"] <= hi_n)]
        if subset.empty:
            continue
        # Representative: node with highest posterior in this size class
        best = subset.loc[subset["posterior"].fillna(0).idxmax()]
        clade_rows.append({
            "clade": label,
            "n_tips": int(best["n_tips"]),
            "height_mean_yrs": best["height_yrs"],
            "hpd_lo_yrs": best["hpd_lo_yrs"],
            "hpd_hi_yrs": best["hpd_hi_yrs"],
            "calendar_mean": best["calendar_mean"],
            "calendar_hpd_lo": best["calendar_hpd_lo"],
            "calendar_hpd_hi": best["calendar_hpd_hi"],
            "posterior": best["posterior"],
        })

    df_clades = pd.DataFrame(clade_rows)
    df_clades.to_csv(RESULTS / "tmrca_by_clade.csv", index=False)

    # Summary text
    lines = ["TMRCA by clade — HPV-16 MCC tree", "=" * 50, f"Reference year: {ref_year}", ""]
    for _, row in df_clades.iterrows():
        cal = row["calendar_mean"]
        lo_c = row["calendar_hpd_lo"]
        hi_c = row["calendar_hpd_hi"]
        lines.append(
            f"{row['clade']}\n"
            f"  Tips          : {row['n_tips']}\n"
            f"  Height (yrs)  : {row['height_mean_yrs']:.1f}  "
            f"[{row['hpd_lo_yrs']:.1f} – {row['hpd_hi_yrs']:.1f}]\n"
            f"  Calendar date : {cal:.0f} AD  [{lo_c:.0f} – {hi_c:.0f} AD]\n"
            f"  Posterior     : {row['posterior']:.3f}\n"
        )
    summary = "\n".join(lines)
    print(summary)
    (RESULTS / "tmrca_summary.txt").write_text(summary, encoding="utf-8")

    # Plot HPD intervals
    fig, ax = plt.subplots(figsize=(9, 5))
    colors = ["#378ADD", "#1D9E75", "#D85A30", "#9B59B6", "#888780"]
    for i, (_, row) in enumerate(df_clades.iterrows()):
        cal = row["calendar_mean"]
        lo_c = row["calendar_hpd_lo"]
        hi_c = row["calendar_hpd_hi"]
        if pd.isna(cal):
            continue
        ax.barh(i, hi_c - lo_c, left=lo_c, height=0.5,
                color=colors[i % len(colors)], alpha=0.6)
        ax.plot(cal, i, "o", color=colors[i % len(colors)], markersize=7)
        ax.text(hi_c + 1, i, f"{cal:.0f} AD", va="center", fontsize=8)

    ax.set_yticks(range(len(df_clades)))
    ax.set_yticklabels(df_clades["clade"], fontsize=9)
    ax.set_xlabel("Calendar year (AD)")
    ax.set_title("TMRCA estimates by clade — 95% HPD intervals")
    ax.axvline(ref_year, color="gray", linestyle="--", linewidth=0.8, label="Most recent sample")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(RESULTS / "tmrca_plot.png", dpi=150)

    print("Saved results/tmrca_nodes.csv")
    print("Saved results/tmrca_by_clade.csv")
    print("Saved results/tmrca_plot.png")
    print("Saved results/tmrca_summary.txt")


if __name__ == "__main__":
    main()
