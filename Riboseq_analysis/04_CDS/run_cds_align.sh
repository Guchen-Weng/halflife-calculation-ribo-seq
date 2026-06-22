#!/bin/bash
# ============================================================
# Script:    run_cds_align.sh
# Pipeline:  Step 04 - CDS Processing
# Purpose:   Align filtered reads to non-overlapping CDS index
#            - Bowtie2 alignment (end-to-end, very-sensitive, -k 1)
#            - Convert SAM to sorted BAM and index
#            - Log alignment statistics
# Input:     ../01-filter/*_clean_final.fastq.gz
#            ../03.1-build-non-overlapping-cds/non_overlapping_cds_index
# Prereq:    conda activate EcoliRNA
# Output:    ../03.2-alignmenttocds/
#            └── *.sorted.bam, *.sorted.bam.bai, logs/
# ============================================================

# Directories
FILTERED_DIR="../01-filter"
REF_DIR="../03.1-build-non-overlapping-cds"
REF_PREFIX="non_overlapping_cds_index"
OUTPUT_DIR="../03.2-alignmenttocds"
LOG_DIR="${OUTPUT_DIR}/logs"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

THREADS=8

# Verify reference index
if [ ! -f "${REF_DIR}/${REF_PREFIX}.1.bt2" ]; then
    echo "Error: CDS index not found in $REF_DIR"
    exit 1
fi

echo "Starting CDS alignment..."
echo "Input directory: $FILTERED_DIR"
echo "Output directory: $OUTPUT_DIR"

# Process files
# Pattern: *_clean_final.fastq.gz
count=0
for r1 in ${FILTERED_DIR}/*_clean_final.fastq.gz; do

    # Extract sample name
    # Filename format: A2_WT_di_1_clean_final.fastq.gz
    # We want: A2_WT_di_1
    basename=$(basename "$r1" _clean_final.fastq.gz)

    echo "--------------------------------------------------"
    echo "Processing $basename..."

    sam_file="${OUTPUT_DIR}/${basename}.sam"
    bam_file="${OUTPUT_DIR}/${basename}.sorted.bam"
    log_file="${LOG_DIR}/${basename}_cds_bowtie2.log"

    # Align to CDS index
    echo "  Aligning to CDS index..."
    bowtie2 -p "$THREADS" -x "${REF_DIR}/${REF_PREFIX}" \
        -U "$r1" \
        --end-to-end --very-sensitive \
        -S "$sam_file" \
        --no-unal \
        -k 1 \
        2> "$log_file"

    if [ $? -ne 0 ]; then
        echo "  Error: Bowtie2 alignment failed for $basename"
        continue
    fi

    # Convert to BAM and sort
    echo "  Converting to sorted BAM..."
    samtools view -bS "$sam_file" | samtools sort -o "$bam_file" -

    if [ $? -ne 0 ]; then
        echo "  Error: Samtools sort failed for $basename"
        rm "$sam_file"
        continue
    fi

    # Index BAM
    echo "  Indexing BAM..."
    samtools index "$bam_file"

    # Cleanup SAM
    rm "$sam_file"

    # Stats
    mapped=$(samtools view -c -F 4 "$bam_file")
    echo "  Mapped reads: $mapped"
    echo "  Mapped reads: $mapped" >> "$log_file"

    ((count++))
done

echo "--------------------------------------------------"
echo "Alignment complete. Processed $count samples."
