#!/usr/bin/env python3
"""
Diagnoses why BA1/BS1 fired incorrectly by replicating InterVar's own
header-matching logic against the real multianno file.

Usage:
    python3 diagnose_ba1.py output/synthetic_htad.intervar.hg38_multianno.txt 48411339
"""
import sys

infile = sys.argv[1]
target_pos = sys.argv[2]

with open(infile) as f:
    lines = f.read().split('\n')

header = lines[0].lstrip('#').split('\t')

Freqs_flgs = {'1000g2015aug_all': 0, 'esp6500siv2_all': 0, 'AF': 0,
              'AF_afr': 0, 'AF_amr': 0, 'AF_eas': 0, 'AF_fin': 0,
              'AF_nfe': 0, 'AF_oth': 0, 'AF_asj': 0}
Allels_flgs = {'Chr': 0, 'Start': 0, 'End': 0, 'Ref': 0, 'Alt': 0}

# Replicate InterVar's exact matching loop (from Intervar.py ~line 1795)
def match_flags(d, header):
    found = {}
    for key in d.keys():
        matched = False
        for i, col in enumerate(header):
            if key == col:
                d[key] = i
                found[key] = i
                matched = True
                break
            if key == "Otherinfo" and (col == key or col == "Otherinfo1"):
                d[key] = i
                found[key] = i
                matched = True
                break
        if not matched:
            found[key] = None  # stayed at default 0 -- NOT actually found
    return found

freq_matches = match_flags(Freqs_flgs, header)
allel_matches = match_flags(Allels_flgs, header)

print("=== Freqs_flgs resolution ===")
for key, idx in freq_matches.items():
    if idx is None:
        print(f"  {key:20s} -> NOT FOUND in header, stays at default index 0 ('{header[0]}')")
    else:
        print(f"  {key:20s} -> index {idx} ('{header[idx]}')")

# Find the target data line
target_line = None
for line in lines[1:]:
    if not line.strip():
        continue
    cls = line.split('\t')
    if len(cls) > 1 and cls[1] == target_pos:
        target_line = cls
        break

if target_line is None:
    print(f"\nNo data line found with Start == {target_pos}")
    sys.exit(1)

print(f"\n=== Values InterVar would actually read for pos {target_pos} ===")
for key, idx in freq_matches.items():
    real_idx = Freqs_flgs[key]
    val = target_line[real_idx] if real_idx < len(target_line) else "<OUT OF RANGE>"
    parseable = "PARSES AS FLOAT" if val.replace('.', '', 1).replace('-', '', 1).isdigit() else "NOT a valid float"
    print(f"  {key:20s} reads cls[{real_idx}] = '{val}'  ({parseable})")

print("\n=== check_BA1 simulation ===")
BA1 = 0
for key in ['1000g2015aug_all', 'esp6500siv2_all', 'AF']:
    idx = Freqs_flgs[key]
    val = target_line[idx] if idx < len(target_line) else None
    try:
        if val is not None and float(val) > 0.05:
            BA1 = 1
            print(f"  -> BA1 SET TO 1 by key '{key}' reading cls[{idx}]='{val}'")
    except ValueError:
        print(f"  -> key '{key}' at cls[{idx}]='{val}' failed float(), no effect (as expected)")
print(f"\nFinal simulated BA1 = {BA1}")
