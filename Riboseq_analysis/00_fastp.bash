#!/bin/bash
# ============================================================
# Script:    00_fastp.bash
# Pipeline:  Step 00 - Raw Read QC
# Purpose:   Quality control of raw paired-end FASTQ files using
#            fastp. Trims adapters, filters low-quality reads
#            (Q20, max 20% unqualified bases, min length 20),
#            and generates HTML/JSON QC reports.
# Input:     ../RawData/<sample>/*_1.fq.gz, *_2.fq.gz
# Output:    ../00-QC/cleaned/<sample>_clean_R{1,2}.fastq.gz
#            ../00-QC/fastp_reports/<sample>_fastp_report.{html,json}
# Author:    Guchen-Weng
# Date:      2025-01
# ============================================================
# Prerequisite: conda activate EcoliRNA

# Project root (one level up from this script's location)
data_dir=".."
raw_dir="../RawData"

# Output directories
output_dir="${data_dir}/00-QC/cleaned"
mkdir -p "$output_dir"

report_dir="${data_dir}/00-QC/fastp_reports"
mkdir -p "$report_dir"

# Log file
LOG_FILE="${data_dir}/00-QC/processing_log_$(date +%Y%m%d_%H%M%S).txt"

# Logging function (writes to both terminal and log file)
log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$LOG_FILE"
}

# Number of threads for fastp
threads=8

# Enter raw data directory
cd "$raw_dir"

log "Starting QC processing in: $raw_dir"

for subdir in */; do
    subdir="${subdir%/}"
    log "Processing subdirectory: $subdir"

    cd "$subdir"

    # Find R1 files in current subdirectory
    for r1_file in *_1.fq.gz; do
        if [ ! -f "$r1_file" ]; then
            log "No R1 files found in $subdir, skipping"
            continue
        fi

        # Derive R2 filename
        r2_file="${r1_file/_1.fq.gz/_2.fq.gz}"

        if [ ! -f "$r2_file" ]; then
            log "No matching R2 file for $r1_file in $subdir, skipping"
            continue
        fi

        # Extract sample base name
        sample_name=$(echo "$r1_file" | sed 's/_1\.fq\.gz//')

        # Output filenames
        output_r1="${output_dir}/${subdir}_clean_R1.fastq.gz"
        output_r2="${output_dir}/${subdir}_clean_R2.fastq.gz"
        output_m="${output_dir}/${subdir}_cleaned_merged.fastq.gz"

        # Report filenames
        html_report="${report_dir}/${subdir}_fastp_report.html"
        json_report="${report_dir}/${subdir}_fastp_report.json"

        log "Processing sample: ${subdir}"
        log "  Input:  $r1_file and $r2_file"
        log "  Output: $output_r1 and $output_r2"

        # Run fastp (paired-end mode)
        fastp -i "$r1_file" -I "$r2_file" \
              -o "$output_r1" -O "$output_r2" \
              -h "$html_report" \
              -j "$json_report" \
              -q 20 -u 20 \
              --detect_adapter_for_pe \
              --n_base_limit 5 \
              --length_required 20 \
              -w "$threads"

        if [ $? -eq 0 ]; then
            log "Successfully completed: ${subdir}"
        else
            log "Processing failed: ${subdir}"
        fi
    done

    cd "$raw_dir"
    log "----------------------------------------"
done

log "All samples processed!"
log "Cleaned reads saved to: $output_dir"
log "QC reports saved to: $report_dir"
