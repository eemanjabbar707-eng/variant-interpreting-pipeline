#!/usr/bin/env python3
"""
08_filter_prioritize.py
Usage: python3 08_filter_prioritize.py <intervar_multianno.txt.intervar> <output.tsv>

Filters and prioritizes variants from an InterVar multianno.intervar file:
  - Deduplicates identical (Chr, Start, Ref, Alt) rows
  - Drops variants classified as Benign / Likely benign
  - Drops variants with gnomAD population AF above a common-variant threshold
  - Ranks remaining variants by ACMG/AMP classification tier, then by
    evidence-code count (more evidence = higher priority within a tier)
"""

import sys
import csv

CLASSIFICATION_RANK = {
    "Pathogenic": 0,
    "Likely pathogenic": 1,
    "Uncertain significance": 2,
    "Likely benign": 3,
    "Benign": 4,
}

GNOMAD_AF_THRESHOLD = 0.01  # 1% - common variant cutoff


def parse_classification(intervar_field):
    """Extract the classification label from the InterVar evidence string."""
    # e.g. "InterVar: Uncertain significance PVS1=0 PS=[...] ..."
    text = intervar_field.replace("InterVar:", "").strip()
    for label in CLASSIFICATION_RANK:
        if text.startswith(label):
            return label
    return "Uncertain significance"


def count_evidence(intervar_field):
    """Rough count of non-zero evidence codes, used as a tiebreaker."""
    count = 0
    for token in ["PVS1", "PS", "PM", "PP", "BA1", "BS", "BP"]:
        idx = intervar_field.find(token + "=")
        if idx == -1:
            continue
        segment = intervar_field[idx: idx + 60]
        count += segment.count("1")
    return count


def safe_float(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return 0.0


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <intervar_multianno.txt.intervar> <output.tsv>")
        sys.exit(1)

    in_path, out_path = sys.argv[1], sys.argv[2]

    print("Filtering and prioritizing variants...")

    with open(in_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)

    intervar_col = next((c for c in rows[0] if "InterVar" in c), None)
    gnomad_col = next((c for c in rows[0] if c.startswith("Freq_gnomAD_genome_ALL")), None)

    if intervar_col is None:
        print("Error: couldn't find an InterVar classification column in the input.")
        sys.exit(1)

    seen = set()
    kept = []
    dropped_dupe = 0
    dropped_benign = 0
    dropped_common = 0

    for row in rows:
        key = (row.get("#Chr") or row.get("Chr"), row["Start"], row["Ref"], row["Alt"])
        if key in seen:
            dropped_dupe += 1
            continue
        seen.add(key)

        classification = parse_classification(row[intervar_col])
        if classification in ("Benign", "Likely benign"):
            dropped_benign += 1
            continue

        af = safe_float(row.get(gnomad_col, ".")) if gnomad_col else 0.0
        if af >= GNOMAD_AF_THRESHOLD:
            dropped_common += 1
            continue

        row["_classification"] = classification
        row["_rank"] = CLASSIFICATION_RANK[classification]
        row["_evidence_count"] = count_evidence(row[intervar_col])
        kept.append(row)

    kept.sort(key=lambda r: (r["_rank"], -r["_evidence_count"]))

    fieldnames = ["#Chr", "Start", "Ref", "Alt", "Ref.Gene", "ExonicFunc.refGene",
                  "_classification", "_evidence_count", gnomad_col or "Freq_gnomAD_genome_ALL",
                  "clinvar: Clinvar"]
    fieldnames = [f for f in fieldnames if f in kept[0] or f.startswith("_")] if kept else fieldnames

    with open(out_path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in kept:
            writer.writerow(row)

    print(f"Total variants read: {len(rows)}")
    print(f"Dropped (duplicate): {dropped_dupe}")
    print(f"Dropped (benign/likely benign): {dropped_benign}")
    print(f"Dropped (common, gnomAD AF >= {GNOMAD_AF_THRESHOLD}): {dropped_common}")
    print(f"Remaining, prioritized: {len(kept)}")
    print(f"Filtering and prioritization completed. Output: {out_path}")


if __name__ == "__main__":
    main()
