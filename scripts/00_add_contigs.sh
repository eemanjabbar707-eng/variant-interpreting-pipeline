#!/bin/bash
# 00_add_contigs.sh
# Adds ##contig header lines to a VCF using the reference FASTA index.
# Usage: bash 00_add_contigs.sh <input.vcf> <output.vcf> <reference.fa>

set -euo pipefail

INPUT_VCF="$1"
OUTPUT_VCF="$2"
REFERENCE="$3"

if [ ! -f "$INPUT_VCF" ]; then
    echo "Error: input file $INPUT_VCF not found" >&2
    exit 1
fi

if [ ! -f "$REFERENCE" ]; then
    echo "Error: reference file $REFERENCE not found" >&2
    exit 1
fi

if [ ! -f "${REFERENCE}.fai" ]; then
    echo "Reference index not found, creating it..."
    samtools faidx "$REFERENCE"
fi

# Build ##contig lines from the .fai index
CONTIG_LINES=$(awk '{printf "##contig=<ID=%s,length=%s>\n", $1, $2}' "${REFERENCE}.fai")

# Insert contig lines right after the ##fileformat line, before any existing
# header lines, avoiding duplicate ##contig entries if some already exist.
awk -v contigs="$CONTIG_LINES" '
    BEGIN { inserted = 0 }
    /^##fileformat/ {
        print
        print contigs
        inserted = 1
        next
    }
    /^##contig=/ { next }  # drop any pre-existing contig lines to avoid duplicates
    { print }
' "$INPUT_VCF" > "$OUTPUT_VCF"

echo "Contig headers added. Output: $OUTPUT_VCF"
