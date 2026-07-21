#!/bin/bash
# run_pipeline.sh
# Usage: ./run_pipeline.sh <input.vcf> <run_name>
#
# Runs the full variant interpretation pipeline end-to-end on a single
# input VCF. Creates all intermediate files under input/<run_name>_stepNN.*
# and a final report at input/<run_name>_final_report.md

set -e  # stop immediately on any step failure

INPUT_VCF=$1
RUN_NAME=$2
REF=reference/GRCh38.fa

if [ -z "$INPUT_VCF" ] || [ -z "$RUN_NAME" ]; then
    echo "Usage: $0 <input.vcf> <run_name>"
    echo "Example: $0 input/synthetic_htad_unannotated.vcf htad"
    exit 1
fi

if [ ! -f "$INPUT_VCF" ]; then
    echo "Error: input file $INPUT_VCF not found"
    exit 1
fi

STEP00="input/${RUN_NAME}_step00.vcf"
STEP02="input/${RUN_NAME}_step02.vcf"
STEP03="input/${RUN_NAME}_step03.vcf"
STEP04_PREFIX="input/${RUN_NAME}_step04"
STEP04="${STEP04_PREFIX}.hg38_multianno.vcf"
STEP05="input/${RUN_NAME}_step05.vcf"
STEP06="input/${RUN_NAME}_step06.vcf"
INTERVAR_PREFIX="input/${RUN_NAME}_intervar_out"
INTERVAR_OUT="${INTERVAR_PREFIX}.hg38_multianno.txt.intervar"
PRIORITIZED="input/${RUN_NAME}_prioritized.tsv"
FINAL_REPORT="input/${RUN_NAME}_final_report.md"

echo "=== [00] Adding contig headers ==="
bash scripts/00_add_contigs.sh "$INPUT_VCF" "$STEP00" "$REF"

echo "=== [01] Validating ==="
bash scripts/01_validate.sh "$STEP00"

echo "=== [02] Normalizing ==="
bash scripts/02_normalize.sh "$STEP00" "$STEP02" "$REF"

echo "=== [03] SnpEff annotation ==="
bash scripts/03_annotate_snpeff.sh "$STEP02" "$STEP03"

echo "=== [04] ClinVar annotation ==="
bash scripts/04_clinvar_annotate.sh "$STEP03" "$STEP04_PREFIX"

echo "=== [05] gnomAD annotation ==="
python3 scripts/05_gnomad_annotate.py "$STEP04" "$STEP05"

echo "=== [06] RepeatMasker annotation ==="
python3 scripts/06_repeatmasker_annotate.py "$STEP05" "$STEP06"

echo "=== [07] InterVar classification (needs 'py2' conda env) ==="
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate py2
bash scripts/07_intervar_classify.sh "$STEP06" "$INTERVAR_PREFIX"
conda deactivate

echo "=== [08] Filtering and prioritizing ==="
python3 scripts/08_filter_prioritize.py "$INTERVAR_OUT" "$PRIORITIZED"

echo "=== [09] Generating report ==="
python3 scripts/09_report.py "$PRIORITIZED" "$FINAL_REPORT"

echo ""
echo "=== Pipeline completed successfully ==="
echo "Final report: $FINAL_REPORT"
