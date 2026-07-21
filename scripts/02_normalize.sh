#!/bin/bash
# 02_normalize.sh
# Usage: ./02_normalize.sh <input.vcf> <output.vcf> <reference.fa>

INPUT=$1
OUTPUT=$2
REF=$3

if [ -z "$INPUT" ] || [ -z "$OUTPUT" ] || [ -z "$REF" ]; then
    echo "Usage: $0 <input.vcf> <output.vcf> <reference.fa>"
    exit 1
fi

if [ ! -f "$INPUT" ]; then
    echo "Error: input file $INPUT not found"
    exit 1
fi

if [ ! -f "$REF" ]; then
    echo "Error: reference file $REF not found"
    exit 1
fi

echo "Normalizing variants with bcftools norm..."

bcftools norm -f "$REF" -c w -Ov -o "$OUTPUT" "$INPUT" 2> normalize.log

MISMATCH_COUNT=$(grep -c "REF_MISMATCH" normalize.log)
if [ "$MISMATCH_COUNT" -gt 0 ]; then
    echo "Warning: $MISMATCH_COUNT REF_MISMATCH record(s) found:"
    grep "REF_MISMATCH" normalize.log
fi

echo "Normalization completed. Output: $OUTPUT"
