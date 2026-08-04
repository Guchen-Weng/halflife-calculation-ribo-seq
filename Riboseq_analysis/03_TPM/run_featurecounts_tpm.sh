#!/bin/bash
# ============================================================
# Script:    run_featurecounts_tpm.sh
# Pipeline:  Step 03 - TPM Calculation
# Purpose:   Run featureCounts on aligned BAM files to generate
#            gene-level counts, then calculate TPM
# Input:     ../02-alignment/*.sorted.bam
# Reference: reference/Escherichia_coli_...gtf
# Output:    03_TPM/counts.txt, 03_TPM/tpm_counts.csv
# Author:    Guchen-Weng
# ============================================================
# Prerequisite: conda activate EcoliRNA

cd "$(dirname "$0")"

GTF="../reference/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.62.chromosome.Chromosome.gtf"

echo "Running featureCounts..."
featureCounts -T 8 \
    -t CDS \
    -g gene_id \
    -a "$GTF" \
    -o counts.txt \
    ../02-alignment/*.sorted.bam

if [ $? -eq 0 ]; then
    echo "featureCounts finished successfully."
    echo "Calculating TPM..."
    python calculate_tpm.py
else
    echo "featureCounts failed!"
    exit 1
fi
