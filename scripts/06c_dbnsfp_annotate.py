#!/usr/bin/env python3
"""
06c_dbnsfp_annotate.py
Usage: python3 06c_dbnsfp_annotate.py <input.vcf> <output.vcf>

Queries the query.genos.us dbNSFP web service per-variant using
DBNSFP_API_EMAIL / DBNSFP_API_KEY from the environment (already saved
in ~/.bashrc), and adds a dbNSFP_score INFO field.
"""

import os
import sys
import time
import requests

DBNSFP_API = "https://query.genos.us/api/dbnsfp"


def fetch_dbnsfp(chrom, pos, ref, alt, email, api_key, max_retries=5):
    chrom_clean = chrom.replace("chr", "")
    params = {
        "email": email,
        "key": api_key,
        "chr": chrom_clean,
        "pos": pos,
        "ref": ref,
        "alt": alt,
    }

    delay = 2.0
    for attempt in range(1, max_retries + 1):
        try:
            resp = requests.get(DBNSFP_API, params=params, timeout=15)
            if resp.status_code == 429:
                print(f"  [rate-limited] {chrom}:{pos} (attempt {attempt}/{max_retries}), waiting {delay:.0f}s...", file=sys.stderr)
                time.sleep(delay)
                delay *= 2
                continue
            resp.raise_for_status()
            data = resp.json()
            score = data.get("CADD_phred") or data.get("score")
            return score
        except Exception as e:
            print(f"  [warn] dbNSFP lookup failed for {chrom}:{pos}: {e}", file=sys.stderr)
            return None

    print(f"  [warn] gave up on {chrom}:{pos} after {max_retries} attempts (still rate-limited)", file=sys.stderr)
    return None


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.vcf> <output.vcf>")
        sys.exit(1)

    in_path, out_path = sys.argv[1], sys.argv[2]

    email = os.environ.get("DBNSFP_API_EMAIL")
    api_key = os.environ.get("DBNSFP_API_KEY")

    if not email or not api_key:
        print("Error: DBNSFP_API_EMAIL / DBNSFP_API_KEY not set in environment.")
        print("Run: source ~/.bashrc")
        sys.exit(1)

    print("Querying dbNSFP (query.genos.us) for pathogenicity scores...")

    with open(in_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            if line.startswith("##"):
                fout.write(line)
                continue
            if line.startswith("#CHROM"):
                fout.write('##INFO=<ID=dbNSFP_score,Number=1,Type=Float,Description="dbNSFP pathogenicity score via query.genos.us">\n')
                fout.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt = fields[0], fields[1], fields[3], fields[4]
            info = fields[7]

            score = fetch_dbnsfp(chrom, pos, ref, alt, email, api_key)
            score_str = score if score is not None else "."
            fields[7] = f"{info};dbNSFP_score={score_str}"

            fout.write("\t".join(fields) + "\n")
            time.sleep(1.0)  # be polite to the API

    print(f"dbNSFP annotation completed. Output: {out_path}")


if __name__ == "__main__":
    main()
