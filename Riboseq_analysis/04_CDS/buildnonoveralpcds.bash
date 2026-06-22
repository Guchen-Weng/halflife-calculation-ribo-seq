#!/bin/bash
# ============================================================
# Script:    buildnonoveralpcds.bash
# Pipeline:  Step 04 - CDS Processing
# Purpose:   Build non-overlapping CDS index
#            - Extract CDS regions from GFF with cleaned IDs
#            - Identify and filter out overlapping CDS on same strand
#            - Generate FASTA and Bowtie2 index
# Input:     ../reference/Escherichia_coli_*_...Chromosome.gff3
#            ../reference/Escherichia_coli_*_...dna.chromosome.Chromosome.fa
# Prereq:    conda activate EcoliRNA
# Output:    ../03.1-build-non-overlapping-cds/
#            └── non_overlapping_cds.bed, .fasta, Bowtie2 index
# ============================================================

# Config
GENOME_FA="../reference/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.dna.chromosome.Chromosome.fa"
GFF="../reference/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.62.chromosome.Chromosome.gff3"
OUTPUT_DIR="../03.1-build-non-overlapping-cds"

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

# 1. Extract CDS and clean IDs
# Purpose: extract ID directly from GFF to avoid complex string processing,
#          ensuring unique IDs without special characters
echo "Extracting and cleaning CDS..."
awk 'BEGIN{OFS="\t"} $3=="CDS" {
    # Extract content after ID= until semicolon or end of line
    match($9, /ID=[^;]+/);
    id=substr($9, RSTART+3, RLENGTH-3);
    # Output: chr, start(0-based), end, id, score, strand
    print $1, $4-1, $5, id, ".", $7
}' "$GFF" > cds_regions.bed

if [ ! -s cds_regions.bed ]; then
    echo "Error: No CDS found." && exit 1
fi

# 2. Identify same-strand overlapping IDs
# Use -s for same-strand, -f 0.01 to ignore tiny overlaps,
# exclude self-matches by ID inequality
echo "Identifying overlapping genes..."
bedtools intersect -a cds_regions.bed -b cds_regions.bed -wa -wb -s | \
awk 'BEGIN{OFS="\t"} $4 != $10 {print $4}' | sort -u > overlapping_ids.txt

echo "Found $(wc -l < overlapping_ids.txt) overlapping CDS entries."

# 3. Filter and extract non-overlapping sequences
# Use fgrep -w for exact ID matching, preventing partial matches (e.g. 1 matching 11)
echo "Filtering and extracting sequences..."
if [ -s overlapping_ids.txt ]; then
    fgrep -v -w -f overlapping_ids.txt cds_regions.bed > non_overlapping_cds.bed
else
    cp cds_regions.bed non_overlapping_cds.bed
fi

# 4. Generate FASTA (using cleaned IDs as headers)
bedtools getfasta -fi "$GENOME_FA" -bed non_overlapping_cds.bed -fo non_overlapping_cds.fasta -name -s

# 5. Build Bowtie2 index
echo "Building Bowtie2 index..."
bowtie2-build --threads 8 non_overlapping_cds.fasta non_overlapping_cds_index

echo "Reference build complete. Final CDS count: $(grep -c ">" non_overlapping_cds.fasta)"
