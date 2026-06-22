#!/bin/bash
# ============================================================
# Script:    build_index.sh
# Purpose:   Build Bowtie2 index from the E. coli K-12 reference
#            genome FASTA included in this directory.
# Usage:     bash build_index.sh
# Output:    Ecoli.1.bt2, Ecoli.2.bt2, Ecoli.3.bt2, Ecoli.4.bt2,
#            Ecoli.rev.1.bt2, Ecoli.rev.2.bt2
# ============================================================

REF_FASTA="Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.dna.chromosome.Chromosome.fa"
INDEX_PREFIX="Ecoli"

if [ ! -f "$REF_FASTA" ]; then
    echo "ERROR: Reference FASTA not found: $REF_FASTA"
    exit 1
fi

echo "Building Bowtie2 index from $REF_FASTA ..."
bowtie2-build "$REF_FASTA" "$INDEX_PREFIX"

if [ $? -eq 0 ]; then
    echo "Index successfully built with prefix: $INDEX_PREFIX"
else
    echo "ERROR: bowtie2-build failed."
    exit 1
fi
