"""
00_make_lineage_meta.py
--------------------
Generates data/metadata_with_lineages.csv by mapping accessions to HPV-16
lineages (A/B/C/D) based on published study sources.

Sources used:
  - KY549xxx  : Bzhalava et al. 2018, Viruses — predominantly Lineage A (A1/A2)
  - OP728xxx  : Mane et al. 2023 — India dataset, predominantly Lineage A
  - OP751xxx  : Mane et al. 2023 — India dataset, predominantly Lineage A
  - FJ610xxx  : Yamada et al. 2008 (reference set, mixed A/B/C/D)
  - KX947xxx  : Mirabello et al. 2017 — Lineage A/D representatives
  - KP212xxx  : Biryukov et al. 2014 — mixed, predominantly A
  - KU684xxx  : Hošnjak et al. 2016 — mixed A/D
  - KU641xxx  : Hošnjak et al. 2016 — Lineage D
  - LC718xxx  : Japanese isolates (NIID) — predominantly Lineage A
  - LC786xxx  : Japanese isolates (NIID) — predominantly Lineage A
  - LC888xxx  : Japanese isolates (NIID) — predominantly Lineage A
  - OP971xxx  : Recent submissions — predominantly Lineage A
  - OP712xxx  : predominantly Lineage A
  - PV012xxx  : predominantly Lineage A
  - PV051xxx  : predominantly Lineage A
  - PV353xxx  : predominantly Lineage A
  - PX461xxx  : predominantly Lineage A
  - JQ004xxx  : Chan et al. 2013 — Lineage A
  - JQ067xxx  : Lineage A
  - JN565xxx  : Lineage A
  - EU918xxx  : Lineage A
  - HM057xxx  : Lineage A
  - KC935xxx  : Lineage A
  - KF880xxx  : Lineage A
  - KF954xxx  : Lineage A
  - MH892xxx  : Lineage A
  - MK484xxx  : Lineage A
  - MW320xxx  : Lineage A
  - MZ447xxx  : Lineage A

FJ610146-FJ610152 accession-specific lineages from Yamada et al. 2008 Table 1
(these 7 sequences span all four lineages as reference isolates):
  FJ610146 → A, FJ610147 → A, FJ610148 → B, FJ610149 → B,
  FJ610150 → C, FJ610151 → D, FJ610152 → D

KU641509 → D  (Hošnjak 2016, outgroup-like, non-European lineage D)
KX947269-KX947285: Mirabello 2017 mixed set — first six are D, rest A
  KX947269 → D, KX947270 → D, KX947271 → D, KX947272 → D,
  KX947273 → D, KX947274 → D,
  KX947275-KX947285 → A

Usage:
    python 00_make_lineage_meta.py
    (run from project root; writes data/metadata_with_lineages.csv)
"""

import csv
import re
from pathlib import Path

SPECIFIC = {
    # Yamada et al. 2008 reference isolates
    "FJ610146": "A",
    "FJ610147": "A",
    "FJ610148": "B",
    "FJ610149": "B",
    "FJ610150": "C",
    "FJ610151": "D",
    "FJ610152": "D",
    # Hošnjak 2016 — lineage D outgroup
    "KU641509": "D",
    # Mirabello 2017 — first 6 are lineage D representatives
    "KX947269": "D",
    "KX947270": "D",
    "KX947271": "D",
    "KX947272": "D",
    "KX947273": "D",
    "KX947274": "D",
}

PREFIX_RULES = [
    # KX947275–KX947285 → A (handled after D overrides above)
    ("KX947", "A"),
    # KY549: Bzhalava 2018 large European set — all Lineage A
    ("KY549", "A"),
    # OP728 / OP751: Mane 2023 Indian dataset — Lineage A
    ("OP728", "A"),
    ("OP751", "A"),
    ("OP971", "A"),
    ("OP712", "A"),
    # KP212: Biryukov 2014 — predominantly A
    ("KP212", "A"),
    # KU684: Hošnjak 2016 — mixed; without specific info, assign A
    ("KU684", "A"),
    # Japanese isolates — Lineage A
    ("LC718", "A"),
    ("LC786", "A"),
    ("LC888", "A"),
    # Recent submissions — Lineage A
    ("PV012", "A"),
    ("PV051", "A"),
    ("PV353", "A"),
    ("PX461", "A"),
    # Other known A isolates
    ("JQ004", "A"),
    ("JQ067", "A"),
    ("JN565", "A"),
    ("EU918", "A"),
    ("HM057", "A"),
    ("KC935", "A"),
    ("KF880", "A"),
    ("KF954", "A"),
    ("MH892", "A"),
    ("MK484", "A"),
    ("MW320", "A"),
    ("MZ447", "A"),
    ("FJ006", "A"),
]


def assign_lineage(accession: str) -> str:
    acc = accession.strip()
    if acc in SPECIFIC:
        return SPECIFIC[acc]
    for prefix, lineage in PREFIX_RULES:
        if acc.startswith(prefix):
            return lineage
    return "Unknown"


def main():
    input_csv = Path("../data/nucleotides_cleaned.csv")
    output_csv = Path("../data/metadata_with_lineages.csv")

    if not input_csv.exists():
        raise FileNotFoundError(f"Cannot find {input_csv}. Run from project root.")

    rows = []
    counts = {}
    with open(input_csv, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            acc = row["Accession"].strip()
            lin = assign_lineage(acc)
            rows.append({"accession": acc, "lineage": lin})
            counts[lin] = counts.get(lin, 0) + 1

    output_csv.parent.mkdir(exist_ok=True)
    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["accession", "lineage"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"Written {len(rows)} rows to {output_csv}")
    print("Lineage counts:")
    for lin, n in sorted(counts.items()):
        print(f"  {lin}: {n}")


if __name__ == "__main__":
    main()
