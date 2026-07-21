#!/bin/bash
# preflight_check.sh
# Usage: bash scripts/preflight_check.sh <input.vcf> <reference.fa>
#
# Catches the exact classes of silent failure we hit manually today:
#   - chromosome naming mismatch vs reference FASTA (e.g. "7" vs "chr7")
#   - missing ##contig header lines
#   - FORMAT/INFO fields used in body but never declared in header
#   - REF allele mismatch against the reference genome at that position
# Exits non-zero with a clear message the moment something is wrong,
# rather than letting a downstream step silently swallow the variant.

set -euo pipefail

VCF="$1"
REF="$2"

fail() { echo "PREFLIGHT FAIL: $1" >&2; exit 1; }
ok()   { echo "PREFLIGHT OK: $1"; }

[ -f "$VCF" ] || fail "Input VCF not found: $VCF"
[ -f "$REF" ] || fail "Reference FASTA not found: $REF"
[ -f "${REF}.fai" ] || samtools faidx "$REF"

N_INPUT=$(grep -vc "^#" "$VCF" || true)
[ "$N_INPUT" -gt 0 ] || fail "Input VCF has zero data rows — nothing to validate."
ok "Found $N_INPUT variant record(s) in input."

# --- 1. Chromosome naming: every CHROM in the VCF must exist in the FASTA index ---
FASTA_CONTIGS=$(cut -f1 "${REF}.fai")
VCF_CHROMS=$(grep -v "^#" "$VCF" | cut -f1 | sort -u)
for c in $VCF_CHROMS; do
  echo "$FASTA_CONTIGS" | grep -qx "$c" || \
    fail "Chromosome '$c' in VCF is not a valid contig name in $REF. Check chr-prefix convention (e.g. '7' vs 'chr7')."
done
ok "All VCF chromosome names match the reference FASTA naming convention."

# --- 2. ##contig header lines must exist for every chromosome used ---
for c in $VCF_CHROMS; do
  grep -q "^##contig=<ID=$c," "$VCF" || \
    fail "Missing '##contig=<ID=$c,...>' header line for chromosome '$c'. bcftools norm will error or silently drop these records."
done
ok "All used chromosomes have proper ##contig header declarations."

# --- 3. Every FORMAT key used in the FORMAT column must have a ##FORMAT= header line ---
FORMAT_KEYS_USED=$(grep -v "^#" "$VCF" | cut -f9 | tr ':' '\n' | sort -u)
for k in $FORMAT_KEYS_USED; do
  grep -q "^##FORMAT=<ID=$k," "$VCF" || \
    fail "FORMAT field '$k' is used in the VCF body but has no '##FORMAT=<ID=$k,...>' header line."
done
ok "All FORMAT fields used are properly declared in the header."

# --- 4. Every INFO key used must have a ##INFO= header line (best-effort parse) ---
INFO_KEYS_USED=$(grep -v "^#" "$VCF" | cut -f8 | tr ';' '\n' | cut -d'=' -f1 | sort -u)
for k in $INFO_KEYS_USED; do
  [ "$k" = "." ] && continue
  grep -q "^##INFO=<ID=$k," "$VCF" || \
    fail "INFO field '$k' is used in the VCF body but has no '##INFO=<ID=$k,...>' header line."
done
ok "All INFO fields used are properly declared in the header."

# --- 5. REF allele sanity check against the actual reference sequence ---
MISMATCHES=0
while IFS=$'\t' read -r chrom pos id ref alt rest; do
  [ "${chrom:0:1}" = "#" ] && continue
  # Only check simple single/multi-base substitutions cleanly (skip symbolic/complex indel edge cases here)
  reflen=${#ref}
  end=$((pos + reflen - 1))
  actual=$(samtools faidx "$REF" "${chrom}:${pos}-${end}" 2>/dev/null | tail -n +2 | tr -d '\n')
  actual_upper=$(echo "$actual" | tr '[:lower:]' '[:upper:]')
  ref_upper=$(echo "$ref" | tr '[:lower:]' '[:upper:]')
  if [ -n "$actual_upper" ] && [ "$actual_upper" != "$ref_upper" ]; then
    echo "  MISMATCH at ${chrom}:${pos} — VCF says REF=$ref, reference genome says $actual" >&2
    MISMATCHES=$((MISMATCHES + 1))
  fi
done < <(grep -v "^#" "$VCF")

[ "$MISMATCHES" -eq 0 ] || fail "$MISMATCHES REF allele mismatch(es) found against $REF (see above). bcftools norm will drop or error on these."
ok "All REF alleles match the reference genome at their stated positions."

# --- 6. Dry-run bcftools norm and confirm variant count is preserved ---
TMP_NORM=$(mktemp --suffix=.vcf)
bcftools norm -f "$REF" "$VCF" -Ov -o "$TMP_NORM" 2>/tmp/preflight_norm.log || \
  fail "bcftools norm itself failed. See /tmp/preflight_norm.log for details."
N_OUTPUT=$(grep -vc "^#" "$TMP_NORM" || true)
rm -f "$TMP_NORM"
[ "$N_OUTPUT" -eq "$N_INPUT" ] || \
  fail "bcftools norm changed variant count: input had $N_INPUT, normalized output has $N_OUTPUT. Investigate before proceeding."
ok "bcftools norm preserves all $N_INPUT variant(s) — safe to proceed."

echo ""
echo "==================================================="
echo " PREFLIGHT PASSED: $VCF is safe to run through the pipeline."
echo "==================================================="
