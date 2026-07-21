#!/usr/bin/env python3
"""
05_gnomad_annotate.py
Usage: python3 05_gnomad_annotate.py <input.vcf> <output.vcf>

Queries the gnomAD GraphQL API per-variant to pull population allele
frequency and adds it as a gnomAD_AF INFO field.
"""

import sys
import time
import requests

GNOMAD_API = "https://gnomad.broadinstitute.org/api"

QUERY = """
query VariantQuery($variantId: String!, $datasetId: DatasetId!) {
  variant(variantId: $variantId, dataset: $datasetId) {
    genome {
      af
    }
    exome {
      af
    }
  }
}
"""

def fetch_gnomad_af(chrom, pos, ref, alt, dataset="gnomad_r4", max_retries=5):
    chrom_clean = chrom.replace("chr", "")
    variant_id = f"{chrom_clean}-{pos}-{ref}-{alt}"
    payload = {
        "query": QUERY,
        "variables": {"variantId": variant_id, "datasetId": dataset},
    }

    delay = 2.0
    for attempt in range(1, max_retries + 1):
        try:
            resp = requests.post(GNOMAD_API, json=payload, timeout=15)
            if resp.status_code == 429:
                print(f"  [rate-limited] {variant_id} (attempt {attempt}/{max_retries}), waiting {delay:.0f}s...", file=sys.stderr)
                time.sleep(delay)
                delay *= 2  # exponential backoff
                continue
            resp.raise_for_status()
            data = resp.json().get("data", {}).get("variant")
            if not data:
                return None
            genome_af = data.get("genome", {}).get("af") if data.get("genome") else None
            exome_af = data.get("exome", {}).get("af") if data.get("exome") else None
            return genome_af if genome_af is not None else exome_af
        except Exception as e:
            print(f"  [warn] gnomAD lookup failed for {variant_id}: {e}", file=sys.stderr)
            return None

    print(f"  [warn] gave up on {variant_id} after {max_retries} attempts (still rate-limited)", file=sys.stderr)
    return None


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.vcf> <output.vcf>")
        sys.exit(1)

    in_path, out_path = sys.argv[1], sys.argv[2]

    print("Querying gnomAD API for population allele frequencies...")

    with open(in_path) as fin, open(out_path, "w") as fout:
        header_written_info = False
        for line in fin:
            if line.startswith("##"):
                fout.write(line)
                continue
            if line.startswith("#CHROM"):
                fout.write('##INFO=<ID=gnomAD_AF,Number=1,Type=Float,Description="gnomAD population allele frequency">\n')
                fout.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt = fields[0], fields[1], fields[3], fields[4]
            info = fields[7]

            af = fetch_gnomad_af(chrom, pos, ref, alt)
            af_str = f"{af:.6g}" if af is not None else "."
            fields[7] = f"{info};gnomAD_AF={af_str}"

            fout.write("\t".join(fields) + "\n")
            time.sleep(1.0)  # be polite to the API

    print(f"gnomAD annotation completed. Output: {out_path}")


if __name__ == "__main__":
    main()
