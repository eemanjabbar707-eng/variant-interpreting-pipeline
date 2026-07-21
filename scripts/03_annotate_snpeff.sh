#!/bin/bash
# 03_annotate_snpeff.sh
# Usage: ./03_annotate_snpeff.sh <input.vcf> <output.vcf>

INPUT=$1
OUTPUT=$2
SNPEFF_JAR=/root/miniconda3/share/snpeff-5.4.0c-0/snpEff.jar
GENOME=GRCh38.99

if [ -z "$INPUT" ] || [ -z "$OUTPUT" ]; then
    echo "Usage: $0 <input.vcf> <output.vcf>"
    exit 1
fi

if [ ! -f "$INPUT" ]; then
    echo "Error: input file $INPUT not found"
    exit 1
fi

echo "Running SnpEff annotation ($GENOME)..."

java -Xmx4g -jar "$SNPEFF_JAR" "$GENOME" "$INPUT" > "$OUTPUT" 2> snpeff.log

if [ $? -ne 0 ]; then
    echo "Error: SnpEff failed. Check snpeff.log"
    exit 1
fi

echo "SnpEff annotation completed. Output: $OUTPUT"
