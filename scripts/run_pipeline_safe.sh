#!/bin/bash
# run_pipeline_safe.sh
# Usage: bash scripts/run_pipeline_safe.sh <input.vcf> <run_name>
#
# Wraps the existing pipeline but checks the variant count after every
# stage. If a stage silently drops variants to zero (or any unexpected
# amount), this aborts immediately with a clear message naming the exact
# stage — instead of letting the pipeline print "completed" on an empty
# file and only failing 2-3 steps later with a cryptic error.

set -euo pipefail

VCF_IN="$1"
RUN="$2"
REF="reference/GRCh38.fa"

count_variants() {
  grep -vc "^#" "$1" 2>/dev/null || echo 0
}

check_step() {
  local label="$1" file="$2" expected="$3"
  [ -f "$file" ] || { echo "ABORT after [$label]: expected output file '$file' does not exist."; exit 1; }
  local n
  n=$(count_variants "$file")
  if [ "$n" -ne "$expected" ]; then
    echo ""
    echo "###########################################################"
    echo "# ABORT after [$label]"
    echo "# Expected $expected variant(s), found $n in: $file"
    echo "# This step silently changed your variant count — stopping"
    echo "# here instead of letting a broken/empty file cascade"
    echo "# through the rest of the pipeline."
    echo "###########################################################"
    exit 1
  fi
  echo "[$label] OK — $n/$expected variant(s) preserved."
}

echo "=== Preflight check ==="
bash scripts/preflight_check.sh "$VCF_IN" "$REF"

N0=$(count_variants "$VCF_IN")

echo ""
echo "=== [00] Adding contig headers ==="
bash scripts/00_add_contigs.sh "$VCF_IN" "input/${RUN}_step00.vcf" "$REF"
check_step "00 add_contigs" "input/${RUN}_step00.vcf" "$N0"

echo ""
echo "=== [01] Validating ==="
bash scripts/01_validate.sh "input/${RUN}_step00.vcf"

echo ""
echo "=== [02] Normalizing ==="
bash scripts/02_normalize.sh "input/${RUN}_step00.vcf" "input/${RUN}_step02.vcf" "$REF"
check_step "02 normalize" "input/${RUN}_step02.vcf" "$N0"

echo ""
echo "=== [03] SnpEff annotation ==="
bash scripts/03_annotate_snpeff.sh "input/${RUN}_step02.vcf" "input/${RUN}_step03.vcf"
check_step "03 snpeff" "input/${RUN}_step03.vcf" "$N0"

echo ""
echo "=== [04] ClinVar annotation ==="
bash scripts/04_clinvar_annotate.sh "input/${RUN}_step03.vcf" "input/${RUN}_step04"
check_step "04 clinvar" "input/${RUN}_step04.hg38_multianno.vcf" "$N0"

echo ""
echo "=== [05] gnomAD annotation ==="
python3 scripts/05_gnomad_annotate.py "input/${RUN}_step04.hg38_multianno.vcf" "input/${RUN}_step05.vcf"
check_step "05 gnomad" "input/${RUN}_step05.vcf" "$N0"

echo ""
echo "=== [06] RepeatMasker annotation ==="
python3 scripts/06_repeatmasker_annotate.py "input/${RUN}_step05.vcf" "input/${RUN}_step06.vcf"
check_step "06 repeatmasker" "input/${RUN}_step06.vcf" "$N0"

echo ""
echo "=== [07] InterVar classification ==="
source activate py2
bash scripts/07_intervar_classify.sh "input/${RUN}_step06.vcf" "input/${RUN}_intervar_out"
conda deactivate
INTERVAR_OUT="input/${RUN}_intervar_out.hg38_multianno.txt.intervar"
[ -f "$INTERVAR_OUT" ] || { echo "ABORT after [07 intervar]: output file missing: $INTERVAR_OUT"; exit 1; }
N_INTERVAR=$(( $(wc -l < "$INTERVAR_OUT") - 1 ))  # minus header row
[ "$N_INTERVAR" -eq "$N0" ] || {
  echo "ABORT after [07 intervar]: expected $N0 classified variant(s), found $N_INTERVAR in $INTERVAR_OUT"
  exit 1
}
echo "[07 intervar] OK — $N_INTERVAR/$N0 variant(s) classified."

echo ""
echo "=== [08] Filtering and prioritizing ==="
python3 scripts/08_filter_prioritize.py "$INTERVAR_OUT" "input/${RUN}_prioritized.tsv"

echo ""
echo "=== [09] Generating report ==="
python3 scripts/09_report.py "input/${RUN}_prioritized.tsv" "input/${RUN}_final_report.md"

echo ""
echo "==================================================="
echo " Pipeline completed successfully — all $N0 variant(s)"
echo " confirmed present at every stage."
echo " Final report: input/${RUN}_final_report.md"
echo "==================================================="
