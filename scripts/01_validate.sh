#!/bin/bash
# 01_validate.sh
# Validates a VCF file: existence, header integrity, sample list, variant count.
# Usage: bash 01_validate.sh <input.vcf>

set -uo pipefail

INPUT_VCF="$1"

echo "========================================="
echo " VariantInterpretationPipeline"
echo " Step 1 : VCF Validation"
echo "========================================="
echo ""
echo "Input:"
echo "$INPUT_VCF"
echo ""

echo "Checking file exists..."
if [ ! -f "$INPUT_VCF" ]; then
    echo "Error: file $INPUT_VCF does not exist" >&2
    exit 1
fi
echo "OK"
echo ""

echo "Checking VCF header..."
# bcftools view -h will fail loudly (and print [W::]/[E::] parse messages)
# if the header or INFO fields are malformed.
if ! bcftools view -h "$INPUT_VCF" > /dev/null 2>/tmp/header_check.log; then
    cat /tmp/header_check.log >&2
    echo "Error: VCF parse error" >&2
    exit 1
fi
echo "Header OK"
echo ""

echo "Checking samples..."
bcftools query -l "$INPUT_VCF"
echo ""

echo "Counting variants..."
VARIANT_COUNT=$(bcftools view -H "$INPUT_VCF" 2>/dev/null | wc -l)
echo "Variants : $VARIANT_COUNT"
echo ""

echo "Generating statistics..."
bcftools stats "$INPUT_VCF" > "${INPUT_VCF%.vcf}.stats.txt" 2>/dev/null
echo ""

echo "Validation completed successfully."
