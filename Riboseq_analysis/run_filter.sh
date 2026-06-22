#!/bin/bash
# ============================================================
# Script:    run_filter.sh
# Pipeline:  Step 01 - tRNA/rRNA Removal
# Purpose:   Remove tRNA and rRNA reads from cleaned FASTQ files
#            using Bowtie2. Performs two sequential alignments:
#            first to tRNA (keeping unaligned), then to rRNA.
# Input:     ../00-QC/cleaned/*_cleaned_merged.fastq.gz
# Reference: reference/tRNA.fa, reference/rRNA.fa
# Output:    ../01-filter/<sample>_clean_final.fastq.gz
# Author:    Guchen-Weng
# Date:      2025-01
# ============================================================
# Prerequisite: conda activate EcoliRNA

# ---- Paths ----
WORKDIR="reference"
READS_DIR="../00-QC/cleaned"
OUTPUT_DIR="../01-filter"

mkdir -p "$OUTPUT_DIR"

# ---- Build Bowtie2 indices if not present ----
cd "$WORKDIR"

if [ ! -f "tRNA.1.bt2" ]; then
    echo "Building tRNA index..."
    bowtie2-build tRNA.fa tRNA
fi

if [ ! -f "rRNA.1.bt2" ]; then
    echo "Building rRNA index..."
    bowtie2-build rRNA.fa rRNA
fi

# ---- Process each sample ----
for r1 in "${READS_DIR}"/*_cleaned_merged.fastq.gz; do
    [ -e "$r1" ] || continue

    basename=$(basename "$r1" _cleaned_merged.fastq.gz)
    echo "Processing $basename..."

    # Step 1: Align to tRNA, keep unaligned reads
    echo "  Aligning to tRNA..."
    bowtie2 --very-sensitive-local -p 8 -x tRNA -U "$r1" \
        --un-gz "${OUTPUT_DIR}/${basename}_no_tRNA.fastq.gz" \
        --no-unal -S /dev/null 2> "${OUTPUT_DIR}/${basename}_tRNA_bowtie2.log"

    # Step 2: Align to rRNA, keep unaligned reads (final clean output)
    echo "  Aligning to rRNA..."
    bowtie2 --very-sensitive-local -p 8 -x rRNA \
        -U "${OUTPUT_DIR}/${basename}_no_tRNA.fastq.gz" \
        --un-gz "${OUTPUT_DIR}/${basename}_clean_final.fastq.gz" \
        --no-unal -S /dev/null 2> "${OUTPUT_DIR}/${basename}_rRNA_bowtie2.log"
done

echo "Done."
