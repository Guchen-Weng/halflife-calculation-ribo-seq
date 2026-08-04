#!/usr/bin/env Rscript
# ============================================================
# Script:    03_read_length_filter.R
# Pipeline:  Step 04 - Monosome (Ribosome Profiling)
# Purpose:   Filter CDS-annotated reads by read length (RFP size
#            selection) to retain monosome-protected fragments.
#            Outputs filtered tables and a per-sample summary of
#            read counts before and after filtering.
# Input:     04_monosome/genome-aligned-bam/<sample>_cds_annotation.csv
# Output:    04_monosome/genome-aligned-bam/<sample>_cds_annotation_len26-32.csv
#            04_monosome/genome-aligned-bam/reads_summary_len26-32.csv
# Parameters: read length 26-32 nt
# Author:    Guchen-Weng
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

# =====================================================================
# 0) Working directory and input files
# =====================================================================
# Run from Riboseq_analysis/ root directory
input_dir  <- "04_monosome/genome-aligned-bam"
output_dir <- "04_monosome/genome-aligned-bam"

stopifnot(dir.exists(input_dir))

files <- c(
  "EV_1_cds_annotation.csv",
  "EV_2_cds_annotation.csv",
  "WT_1_cds_annotation.csv",
  "WT_2_cds_annotation.csv"
)

# Build full paths and verify existence
file_paths <- file.path(input_dir, files)
stopifnot(all(file.exists(file_paths)))

# =====================================================================
# 1) Read-length filtering range
# =====================================================================
len_min <- 26L
len_max <- 32L

# =====================================================================
# 2) Batch filtering and read-count summary
# =====================================================================
stats <- rbindlist(lapply(file_paths, function(fp) {
  f <- basename(fp)
  dt <- fread(fp)

  # Convert relevant columns to numeric types
  dt[, read_length := as.integer(read_length)]
  dt[, readcount   := as.numeric(readcount)]

  # Total read counts before length filtering
  reads_before <- dt[, sum(readcount, na.rm = TRUE)]

  # Retain reads within the monosome-protected range
  dt_filt <- dt[read_length >= len_min & read_length <= len_max]
  reads_after <- dt_filt[, sum(readcount, na.rm = TRUE)]

  # Write filtered table
  out_file <- file.path(
    output_dir,
    str_replace(f, "\\.csv$", paste0("_len", len_min, "-", len_max, ".csv"))
  )
  fwrite(dt_filt, out_file)

  # Return per-sample filtering statistics
  data.table(
    sample       = str_replace(f, "_cds_annotation\\.csv$", ""),
    in_file      = f,
    out_file     = basename(out_file),
    reads_before = reads_before,
    reads_after  = reads_after,
    frac_kept    = fifelse(reads_before > 0, reads_after / reads_before, NA_real_)
  )
}), fill = TRUE)

# =====================================================================
# 3) Save summary table
# =====================================================================
summary_file <- file.path(output_dir, paste0("reads_summary_len", len_min, "-", len_max, ".csv"))
fwrite(stats, summary_file)

message("Filtering summary:")
print(stats)

message("Done. Filtered CSV files and summary table saved in: ", output_dir)
message("Summary file: ", basename(summary_file))
