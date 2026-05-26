"""
05_branch_outliers.py
---------------------
Detects long-branch outliers in the HPV-16 MCC tree.

Long branches can indicate:
  - Misaligned or chimeric sequences
  - Genuinely divergent lineages (e.g. HPV-16 B/C/D vs A)
  - Sequences with excess ambiguous bases that inflated divergence

Flags tips and internal nodes whose branch length exceeds:
  - 2 SD above the mean (warning)
  - 3 SD above the mean (severe — review strongly recommended)

Outputs:
  - results/branch_lengths.csv           (all branches with z-scores)
  - results/branch_outliers.csv          (flagged branches only)
  - results/branch_length_dist.png       (distribution + thresholds)
  - results/branch_outlier_summary.txt   (human-readable report)

Usage:
    python code/05_branch_outliers.py --tree code/<your>.tree
"""

import argparse
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
    p.add_argument("--warn_sd", type=float, default=2.0,
                   help="SD threshold for warning. Default: 2.0")
    p.add_argument("--severe_sd", type=float, default=3.0,
                   help="SD threshold for severe outlier. Default: 3.0")
    return p.parse_args()


def collect_branches(tree) -> list[dict]:
    rows = []
    for clade in tree.find_clades():
        bl = clade.branch_length
        if bl is None:
            continue
        is_tip = clade.is_terminal()
        name = clade.name if is_tip else f"internal_{id(clade)}"
        n_descendants = len(clade.get_terminals())
        rows.append({
            "name": name,
            "branch_length": bl,
            "is_terminal": is_tip,
            "n_descendants": n_descendants,
        })
    return rows


def main():
    args = parse_args()
    tree_path = Path(args.tree)
    if not tree_path.exists():
        sys.exit(f"Tree not found: {tree_path}")

    print(f"Loading tree: {tree_path}")
    tree = Phylo.read(str(tree_path), "nexus")

    rows = collect_branches(tree)
    df = pd.DataFrame(rows)
    print(f"  {len(df)} branches collected "
          f"({df['is_terminal'].sum()} terminal, {(~df['is_terminal']).sum()} internal)")

    # Compute z-scores separately for terminal and internal branches
    for group, mask in [("terminal", df["is_terminal"]), ("internal", ~df["is_terminal"])]:
        sub = df.loc[mask, "branch_length"]
        mean, std = sub.mean(), sub.std()
        df.loc[mask, f"{group}_mean"] = mean
        df.loc[mask, f"{group}_std"] = std
        df.loc[mask, "z_score"] = (sub - mean) / std

    # Flag thresholds
    df["flag"] = "ok"
    df.loc[df["z_score"] > args.warn_sd, "flag"] = "warning"
    df.loc[df["z_score"] > args.severe_sd, "flag"] = "severe"

    df.sort_values("z_score", ascending=False, inplace=True)
    df.to_csv(RESULTS / "branch_lengths.csv", index=False)

    outliers = df[df["flag"] != "ok"]
    outliers.to_csv(RESULTS / "branch_outliers.csv", index=False)

    warn = df[df["flag"] == "warning"]
    severe = df[df["flag"] == "severe"]
    terminals = df[df["is_terminal"]]
    t_mean = terminals["branch_length"].mean()
    t_std = terminals["branch_length"].std()
    t_median = terminals["branch_length"].median()

    summary_lines = [
        "Branch length outlier report — HPV-16 MCC tree",
        "=" * 55,
        f"Total branches     : {len(df)}",
        f"  Terminal         : {df['is_terminal'].sum()}",
        f"  Internal         : {(~df['is_terminal']).sum()}",
        "",
        "Terminal branch statistics:",
        f"  Mean             : {t_mean:.4f} yrs",
        f"  Std dev          : {t_std:.4f} yrs",
        f"  Median           : {t_median:.4f} yrs",
        f"  Warning (>{args.warn_sd}σ)  : {len(warn)} branches",
        f"  Severe (>{args.severe_sd}σ)   : {len(severe)} branches",
        "",
        "Severe outliers (review these sequences):",
    ]
    for _, row in severe.iterrows():
        summary_lines.append(
            f"  {row['name']:<40} BL={row['branch_length']:.4f}  z={row['z_score']:.2f}"
            + ("  [TIP]" if row["is_terminal"] else "  [internal]")
        )

    if severe.empty:
        summary_lines.append("  None — no severe outliers detected.")

    summary = "\n".join(summary_lines)
    print(summary)
    (RESULTS / "branch_outlier_summary.txt").write_text(summary)

    # Plot: distribution + thresholds
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for ax, (label, mask) in zip(
        axes, [("Terminal branches", df["is_terminal"]), ("Internal branches", ~df["is_terminal"])]
    ):
        sub = df.loc[mask, "branch_length"]
        mean, std = sub.mean(), sub.std()
        ax.hist(sub, bins=50, color="#378ADD", alpha=0.7, edgecolor="white", linewidth=0.3)
        ax.axvline(mean, color="#333", linewidth=1.2, linestyle="-", label=f"Mean={mean:.2f}")
        ax.axvline(mean + args.warn_sd * std, color="#EF9F27", linewidth=1.2,
                   linestyle="--", label=f"+{args.warn_sd}σ (warning)")
        ax.axvline(mean + args.severe_sd * std, color="#D85A30", linewidth=1.2,
                   linestyle="--", label=f"+{args.severe_sd}σ (severe)")
        ax.set_xlabel("Branch length (years)")
        ax.set_ylabel("Count")
        ax.set_title(f"{label} — HPV-16 MCC tree")
        ax.legend(fontsize=8)

    fig.suptitle("Branch length distribution and outlier thresholds", fontsize=11)
    fig.tight_layout()
    fig.savefig(RESULTS / "branch_length_dist.png", dpi=150)

    # Log-scale version
    fig2, ax2 = plt.subplots(figsize=(8, 5))
    all_bl = df["branch_length"].values
    log_bl = np.log10(all_bl[all_bl > 0])
    ax2.hist(log_bl, bins=60, color="#1D9E75", alpha=0.7, edgecolor="white", linewidth=0.3)
    ax2.set_xlabel("log₁₀(branch length)")
    ax2.set_ylabel("Count")
    ax2.set_title("Branch length distribution — log scale (reveals multimodality)")
    fig2.tight_layout()
    fig2.savefig(RESULTS / "branch_length_dist_log.png", dpi=150)

    print(f"\nOutliers: {len(warn)} warnings, {len(severe)} severe")
    print("Saved results/branch_lengths.csv")
    print("Saved results/branch_outliers.csv")
    print("Saved results/branch_length_dist.png")
    print("Saved results/branch_length_dist_log.png")
    print("Saved results/branch_outlier_summary.txt")


if __name__ == "__main__":
    main()
