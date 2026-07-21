#!/bin/bash
# 07_intervar_classify.sh
# Usage: ./07_intervar_classify.sh <input.vcf> <output_prefix>

INPUT=$1
OUT_PREFIX=$2
INTERVAR_DIR=software/InterVar
ANNOVAR_DIR=software/annovar
HUMANDB=$ANNOVAR_DIR/humandb

if [ -z "$INPUT" ] || [ -z "$OUT_PREFIX" ]; then
    echo "Usage: $0 <input.vcf> <output_prefix>"
    exit 1
fi

if [ ! -f "$INPUT" ]; then
    echo "Error: input file $INPUT not found"
    exit 1
fi

echo "Running InterVar ACMG/AMP classification..."

python2 "$INTERVAR_DIR/Intervar.py" \
    -b hg38 \
    -i "$INPUT" \
    --input_type=VCF \
    -o "$OUT_PREFIX" \
    -d "$HUMANDB" \
    --table_annovar="$ANNOVAR_DIR/table_annovar.pl" \
    --convert2annovar="$ANNOVAR_DIR/convert2annovar.pl" \
    --annotate_variation="$ANNOVAR_DIR/annotate_variation.pl" \
    -c "$INTERVAR_DIR/config.ini" 2> intervar.log

if [ $? -ne 0 ]; then
    echo "Error: InterVar failed. Check intervar.log"
    exit 1
fi

echo "InterVar classification completed. Output: ${OUT_PREFIX}.hg38_multianno.txt.intervar"
