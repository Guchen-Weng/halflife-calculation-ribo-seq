#!/bin/bash
# ============================================================
# Script:    02_alignment.sh
# Pipeline:  Step 02 - Bowtie2 Genome Alignment
# Purpose:   Align filtered Ribo-seq reads to E. coli K-12 genome
#            using Bowtie2. Converts SAM to sorted/indexed BAM.
# Input:     ../01-filter/*_clean_final.fastq.gz
# Reference: reference/Ecoli (Bowtie2 index)
# Output:    ../02-alignment/<sample>.sorted.bam (.bai)
# Author:    Guchen-Weng
# ============================================================
# Prerequisite: conda activate EcoliRNA

# ---- Paths ----
FILTERED_DIR="../01-filter"
REF_DIR="reference"
REF_PREFIX="Ecoli"
OUTPUT_DIR="../02-alignment"
LOG_DIR="${OUTPUT_DIR}/logs"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

THREADS=8

# ---- Validate reference index ----
if [ ! -f "${REF_DIR}/${REF_PREFIX}.1.bt2" ]; then
    echo "ERROR: Bowtie2 index not found in $REF_DIR"
    echo "Build it with: cd reference && bash build_index.sh"
    exit 1
fi

echo "Starting genome alignment..."
echo "Input directory: $FILTERED_DIR"
echo "Output directory: $OUTPUT_DIR"

# ---- Process each sample ----
count=0
for r1 in ${FILTERED_DIR}/*_clean_final.fastq.gz; do

    basename=$(basename "$r1" _clean_final.fastq.gz)

    echo "--------------------------------------------------"
    echo "Processing $basename..."

    sam_file="${OUTPUT_DIR}/${basename}.sam"
    bam_file="${OUTPUT_DIR}/${basename}.sorted.bam"
    log_file="${LOG_DIR}/${basename}_genome_bowtie2.log"

    # Bowtie2 alignment (end-to-end, keep best alignment only)
    echo "  Aligning to E. coli genome..."
    bowtie2 -p "$THREADS" -x "${REF_DIR}/${REF_PREFIX}" \
        -U "$r1" \
        --end-to-end --very-sensitive \
        -S "$sam_file" \
        --no-unal \
        -k 1 \
        2> "$log_file"

    if [ $? -ne 0 ]; then
        echo "  ERROR: Bowtie2 alignment failed for $basename"
        continue
    fi

    # Convert SAM to sorted BAM
    echo "  Converting to sorted BAM..."
    samtools view -bS "$sam_file" | samtools sort -o "$bam_file" -

    if [ $? -ne 0 ]; then
        echo "  ERROR: Samtools sort failed for $basename"
        rm "$sam_file"
        continue
    fi

    # Index BAM
    echo "  Indexing BAM..."
    samtools index "$bam_file"

    # Remove intermediate SAM
    rm "$sam_file"

    # Alignment statistics
    mapped=$(samtools view -c -F 4 "$bam_file")
    echo "  Mapped reads: $mapped"
    echo "  Mapped reads: $mapped" >> "$log_file"

    ((count++))
done

echo "--------------------------------------------------"
echo "Alignment complete. Processed $count samples."
