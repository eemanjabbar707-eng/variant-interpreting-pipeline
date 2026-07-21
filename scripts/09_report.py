#!/usr/bin/env python3
"""
09_report.py
Usage: python3 09_report.py <prioritized_variants.tsv> <output.md>

Generates a Markdown summary report from the filtered/prioritized
variant table: overview counts, a per-classification breakdown, and
a full variant table sorted by priority.
"""

import sys
import csv
from collections import Counter
from datetime import datetime


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <prioritized_variants.tsv> <output.md>")
        sys.exit(1)

    in_path, out_path = sys.argv[1], sys.argv[2]

    print("Generating variant interpretation report...")

    with open(in_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)

    class_counts = Counter(r["_classification"] for r in rows)
    gene_counts = Counter(r.get("Ref.Gene", ".") for r in rows)

    lines = []
    lines.append("# Variant Interpretation Report")
    lines.append("")
    lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    lines.append("")
    lines.append("## Overview")
    lines.append("")
    lines.append(f"- Total variants (post-filter): **{len(rows)}**")
    lines.append(f"- Genes represented: **{len(gene_counts)}** ({', '.join(sorted(gene_counts))})")
    lines.append("")
    lines.append("## Classification breakdown")
    lines.append("")
    lines.append("| Classification | Count |")
    lines.append("|---|---|")
    for label in ["Pathogenic", "Likely pathogenic", "Uncertain significance",
                  "Likely benign", "Benign"]:
        if class_counts.get(label):
            lines.append(f"| {label} | {class_counts[label]} |")
    lines.append("")

    priority = [r for r in rows if r["_classification"] in ("Pathogenic", "Likely pathogenic")]
    if priority:
        lines.append("## Variants of interest (Pathogenic / Likely pathogenic)")
        lines.append("")
        lines.append("| Chr | Pos | Ref | Alt | Gene | Classification | Evidence count |")
        lines.append("|---|---|---|---|---|---|---|")
        for r in priority:
            lines.append(
                f"| {r.get('#Chr','.')} | {r.get('Start','.')} | {r.get('Ref','.')} | "
                f"{r.get('Alt','.')} | {r.get('Ref.Gene','.')} | {r['_classification']} | "
                f"{r.get('_evidence_count','.')} |"
            )
        lines.append("")
    else:
        lines.append("## Variants of interest (Pathogenic / Likely pathogenic)")
        lines.append("")
        lines.append("None found in this run.")
        lines.append("")

    lines.append("## Full variant table (prioritized)")
    lines.append("")
    lines.append("| Chr | Pos | Ref | Alt | Gene | Classification | Evidence count | gnomAD AF |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for r in rows:
        gnomad_col = next((k for k in r if k.startswith("Freq_gnomAD_genome_ALL")), None)
        af = r.get(gnomad_col, ".") if gnomad_col else "."
        lines.append(
            f"| {r.get('#Chr','.')} | {r.get('Start','.')} | {r.get('Ref','.')} | "
            f"{r.get('Alt','.')} | {r.get('Ref.Gene','.')} | {r['_classification']} | "
            f"{r.get('_evidence_count','.')} | {af} |"
        )
    lines.append("")

    with open(out_path, "w") as out:
        out.write("\n".join(lines))

    print(f"Report generation completed. Output: {out_path}")


if __name__ == "__main__":
    main()
