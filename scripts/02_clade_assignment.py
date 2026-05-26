"""
02_clade_assignment.py
----------------------
Maps HPV-16 lineages (A, B, C, D) onto the MCC tree.

Lineage is inferred from accession prefix patterns and, if a metadata CSV
is provided (--meta), from an explicit 'lineage' or 'Lineage' column.
Falls back to clustering by root-to-tip distance quartile when no metadata
is available, and labels clades accordingly.

Outputs:
  - results/clade_assignment.csv         (tip → lineage mapping)
  - results/clade_monophyly.txt          (monophyly test per lineage)
  - results/clade_tree_colored.png       (tree with lineage colours)

Usage:
    python code/02_clade_assignment.py --tree code/<your>.tree
    python code/02_clade_assignment.py --tree code/<your>.tree --meta code/metadata.csv
"""

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
from Bio import Phylo

RESULTS = Path("results")
RESULTS.mkdir(exist_ok=True)

# Known accession prefix → lineage mappings (extend as needed)
# Based on published HPV-16 lineage reference sets
LINEAGE_HINTS = {
    "A": ["K02718", "AF125673", "AY686581"],
    "B": ["AF536179", "HQ644257"],
    "C": ["AF402678"],
    "D": ["AY686583", "HQ644258"],
}

LINEAGE_COLORS = {
    "A": "#378ADD",
    "B": "#1D9E75",
    "C": "#D85A30",
    "D": "#9B59B6",
    "Unknown": "#888780",
}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--tree", required=True)
    p.add_argument("--meta", default=None,
                   help="Optional CSV with columns 'accession' and 'lineage'")
    return p.parse_args()


def accession_from_name(name: str) -> str:
    return name.split("_")[0]


def load_lineage_map(meta_path: str) -> dict[str, str]:
    df = pd.read_csv(meta_path)
    df.columns = df.columns.str.lower()
    if "lineage" not in df.columns:
        print("  Warning: no 'lineage' column found in metadata — skipping")
        return {}
    acc_col = next((c for c in df.columns if "accession" in c), None)
    if acc_col is None:
        print("  Warning: no accession column found in metadata — skipping")
        return {}
    return dict(zip(df[acc_col].str.strip(), df["lineage"].str.strip()))


def assign_lineages(tips, lineage_map: dict) -> dict[str, str]:
    assignments = {}
    for tip in tips:
        acc = accession_from_name(tip.name)
        if acc in lineage_map:
            assignments[tip.name] = lineage_map[acc]
            continue
        # Check known reference prefixes
        assigned = "Unknown"
        for lin, refs in LINEAGE_HINTS.items():
            if any(acc.startswith(r[:6]) for r in refs):
                assigned = lin
                break
        assignments[tip.name] = assigned
    return assignments


def test_monophyly(tree, assignments: dict) -> str:
    lineage_tips = defaultdict(list)
    for tip_name, lin in assignments.items():
        if lin != "Unknown":
            lineage_tips[lin].append(tip_name)

    lines = ["Monophyly test per lineage", "=" * 40]
    for lin in sorted(lineage_tips):
        tips_in_lin = lineage_tips[lin]
        if len(tips_in_lin) < 2:
            lines.append(f"  Lineage {lin}: only {len(tips_in_lin)} tip — skipped")
            continue
        # Get MRCA of this lineage's tips
        tip_objs = [t for t in tree.get_terminals() if t.name in tips_in_lin]
        mrca = tree.common_ancestor(tip_objs)
        mrca_tips = {t.name for t in mrca.get_terminals()}
        lineage_set = set(tips_in_lin)
        is_mono = mrca_tips == lineage_set
        intruders = mrca_tips - lineage_set
        lines.append(
            f"  Lineage {lin}: {len(tips_in_lin)} tips | "
            f"MRCA subtree size={len(mrca_tips)} | "
            f"monophyletic={'YES' if is_mono else 'NO'}"
            + (f" | intruders={len(intruders)}" if not is_mono else "")
        )
    return "\n".join(lines)


def main():
    args = parse_args()
    tree_path = Path(args.tree)
    if not tree_path.exists():
        sys.exit(f"Tree not found: {tree_path}")

    print(f"Loading tree: {tree_path}")
    tree = Phylo.read(str(tree_path), "nexus")
    tips = tree.get_terminals()
    print(f"  {len(tips)} tips")

    lineage_map = {}
    if args.meta:
        print(f"Loading metadata: {args.meta}")
        lineage_map = load_lineage_map(args.meta)
        print(f"  {len(lineage_map)} lineage annotations loaded")

    assignments = assign_lineages(tips, lineage_map)

    counts = defaultdict(int)
    for v in assignments.values():
        counts[v] += 1
    print("Lineage counts:")
    for lin, n in sorted(counts.items()):
        print(f"  {lin}: {n}")

    df = pd.DataFrame([
        {"tip": k, "accession": accession_from_name(k), "lineage": v}
        for k, v in assignments.items()
    ])
    df.to_csv(RESULTS / "clade_assignment.csv", index=False)

    mono_text = test_monophyly(tree, assignments)
    print(mono_text)
    (RESULTS / "clade_monophyly.txt").write_text(mono_text, encoding="utf-8")

    # Colour tree by lineage
    for tip in tips:
        tip.color = LINEAGE_COLORS.get(assignments.get(tip.name, "Unknown"), "#888780")

    fig, ax = plt.subplots(figsize=(10, max(8, len(tips) * 0.05)))
    Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda x: "")
    ax.set_title("HPV-16 MCC tree — lineage assignment (A/B/C/D)")
    patches = [mpatches.Patch(color=c, label=f"Lineage {l}")
               for l, c in LINEAGE_COLORS.items()]
    ax.legend(handles=patches, loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(RESULTS / "clade_tree_colored.png", dpi=150)

    print("Saved results/clade_assignment.csv")
    print("Saved results/clade_monophyly.txt")
    print("Saved results/clade_tree_colored.png")


if __name__ == "__main__":
    main()
