#!/bin/bash
# ============================================================
# Script:    01-alignment-double.bash
# Pipeline:  Step 01 - Read Alignment
# Purpose:   Align cleaned paired-end FASTQ reads to the E. coli
#            K-12 reference genome using Bowtie2. Converts SAM to
#            sorted/indexed BAM and computes mapping statistics.
# Input:     ../00-QC/cleaned/*_R{1,2}.fastq.gz
#            reference/Ecoli (Bowtie2 index; build with reference/build_index.sh)
# Output:    ../01-alignment/<sample>.sorted.bam (.bai)
#            ../01-alignment/logs/<sample>_bowtie2.log
# Author:    Guchen-Weng
# ============================================================

# ---- Paths ----
clean_fastq_dir="../00-QC/cleaned"
ref_dir="./reference"  # Reference genome directory (see reference/README.md)

# Bowtie2 index prefix (without extension)
ref_prefix="Ecoli"

output_dir="../01-alignment"
mkdir -p "$output_dir"

# Log directory
log_dir="${output_dir}/logs"
mkdir -p "$log_dir"

threads=8

# ---- Validate reference index ----
if [ ! -f "${ref_dir}/${ref_prefix}.1.bt2" ]; then
    echo "ERROR: Bowtie2 index file ${ref_prefix}.1.bt2 not found in ${ref_dir}"
    echo "Please ensure the index is built and ref_prefix is set correctly."
    exit 1
fi

# ---- Main alignment loop ----
cd "$clean_fastq_dir"

total_pairs=0
processed_pairs=0

echo "Starting paired-end FASTQ processing..."

r1_files=(*_R1.fastq.gz)
total_pairs=${#r1_files[@]}

if [ $total_pairs -eq 0 ]; then
    echo "ERROR: No R1 files found in $clean_fastq_dir!"
    exit 1
fi

echo "Found $total_pairs FASTQ pairs to process."

for r1_file in "${r1_files[@]}"; do
    # Derive R2 filename
    r2_file="${r1_file/_clean_R1.fastq.gz/_clean_R2.fastq.gz}"

    if [ ! -f "$r2_file" ]; then
        echo "ERROR: No matching R2 file: $r2_file, skipping this sample"
        continue
    fi

    # Extract sample base name
    sample_name=$(echo "$r1_file" | sed 's/_R1\.fastq\.gz//')

    # Output filenames
    sam_file="${output_dir}/${sample_name}.sam"
    bam_file="${output_dir}/${sample_name}.sorted.bam"
    log_file="${log_dir}/${sample_name}_bowtie2.log"

    echo "Processing sample: $sample_name"
    echo "  Input:  $r1_file and $r2_file"
    echo "  Output: $bam_file"

    # Step 1: Bowtie2 alignment (paired-end)
    echo "  Running Bowtie2 alignment..."
    bowtie2 -x "${ref_dir}/${ref_prefix}" \
            -1 "$r1_file" \
            -2 "$r2_file" \
            -S "$sam_file" \
            --no-unal \
            --un-conc "${output_dir}/${sample_name}_unaligned.fastq.gz" \
            --threads "$threads" > "$log_file" 2>&1

    if [ $? -ne 0 ]; then
        echo "ERROR: Bowtie2 alignment failed for $sample_name! Check log: $log_file"
        continue
    fi

    # Step 2: SAM -> sorted BAM
    echo "  Converting SAM to sorted BAM..."
    samtools view -bS "$sam_file" | samtools sort -o "$bam_file" -

    if [ $? -ne 0 ]; then
        echo "ERROR: SAM-to-BAM conversion failed for $sample_name!"
        continue
    fi

    # Step 3: Index BAM
    echo "  Creating BAM index..."
    samtools index "$bam_file"

    # Step 4: Remove intermediate SAM to save space
    rm "$sam_file"

    # Step 5: Compute mapping rate
    echo "  Computing mapping statistics..."
    total_reads=$(samtools view -c "$bam_file")
    mapped_reads=$(samtools view -c -F 4 "$bam_file")
    if [ $total_reads -gt 0 ]; then
        mapping_rate=$(echo "scale=2; $mapped_reads * 100 / $total_reads" | bc)
        echo "  Mapping rate: $mapping_rate% ($mapped_reads/$total_reads)"
        echo "Mapping stats: $mapping_rate% ($mapped_reads/$total_reads)" >> "$log_file"
    else
        echo "  WARNING: No reads found!"
        echo "WARNING: No reads found!" >> "$log_file"
    fi

    ((processed_pairs++))
    echo "  Successfully completed: $sample_name"
    echo "----------------------------------------"
done

# ---- Final report ----
echo "All samples processed!"
echo "Successfully completed: $processed_pairs/$total_pairs pairs"
echo "Alignment results saved to: $output_dir"
echo "Log files saved to: $log_dir"
