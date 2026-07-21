#!/usr/bin/env python3
"""
06_repeatmasker_annotate.py
Usage: python3 06_repeatmasker_annotate.py <input.vcf> <output.vcf>

Queries the UCSC REST API's rmsk (RepeatMasker) track per-variant to flag
whether a position falls inside a repetitive/low-complexity region.
"""

import sys
import time
import requests

UCSC_API = "https://api.genome.ucsc.edu/getData/track"
GENOME = "hg38"
TRACK = "rmsk"


def fetch_repeat_flag(chrom, pos, max_retries=5):
    params = {
        "genome": GENOME,
        "track": TRACK,
        "chrom": chrom,
        "start": int(pos) - 1,  # UCSC is 0-based
        "end": int(pos),
    }

    delay = 2.0
    for attempt in range(1, max_retries + 1):
        try:
            resp = requests.get(UCSC_API, params=params, timeout=15)
            if resp.status_code == 429:
                print(f"  [rate-limited] {chrom}:{pos} (attempt {attempt}/{max_retries}), waiting {delay:.0f}s...", file=sys.stderr)
                time.sleep(delay)
                delay *= 2
                continue
            resp.raise_for_status()
            data = resp.json()
            track_data = data.get(TRACK, [])
            if track_data:
                repeat_name = track_data[0].get("repName", "unknown")
                repeat_class = track_data[0].get("repClass", "unknown")
                return f"{repeat_class}:{repeat_name}"
            return None
        except Exception as e:
            print(f"  [warn] RepeatMasker lookup failed for {chrom}:{pos}: {e}", file=sys.stderr)
            return None

    print(f"  [warn] gave up on {chrom}:{pos} after {max_retries} attempts (still rate-limited)", file=sys.stderr)
    return None


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.vcf> <output.vcf>")
        sys.exit(1)

    in_path, out_path = sys.argv[1], sys.argv[2]

    print("Querying UCSC RepeatMasker track for repeat-region overlap...")

    with open(in_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            if line.startswith("##"):
                fout.write(line)
                continue
            if line.startswith("#CHROM"):
                fout.write('##INFO=<ID=RepeatMasker,Number=1,Type=String,Description="Overlapping repeat class:name from UCSC rmsk track, if any">\n')
                fout.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            chrom, pos = fields[0], fields[1]
            info = fields[7]

            repeat_flag = fetch_repeat_flag(chrom, pos)
            flag_str = repeat_flag if repeat_flag is not None else "."
            fields[7] = f"{info};RepeatMasker={flag_str}"

            fout.write("\t".join(fields) + "\n")
            time.sleep(1.0)  # be polite to the API

    print(f"RepeatMasker annotation completed. Output: {out_path}")


if __name__ == "__main__":
    main()
