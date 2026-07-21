#!/bin/bash
# 04_clinvar_annotate.sh
# Usage: ./04_clinvar_annotate.sh <input.vcf> <output_prefix>

INPUT=$1
OUT_PREFIX=$2
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

echo "Annotating with ClinVar (clinvar_20240917)..."

perl "$ANNOVAR_DIR/table_annovar.pl" "$INPUT" "$HUMANDB" \
    -buildver hg38 \
    -out "$OUT_PREFIX" \
    -remove \
    -protocol clinvar_20240917 \
    -operation f \
    -nastring . \
    -vcfinput 2> clinvar.log

if [ $? -ne 0 ]; then
    echo "Error: ClinVar annotation failed. Check clinvar.log"
    exit 1
fi

echo "ClinVar annotation completed. Output: ${OUT_PREFIX}.hg38_multianno.vcf"
