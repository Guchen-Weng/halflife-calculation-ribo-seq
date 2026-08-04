#!/usr/bin/env Rscript
# ============================================================
# Script:    02_p_site_calculation.R
# Pipeline:  Step 04 - Monosome (Ribosome Profiling)
# Purpose:   Assign P/A/E-site genomic coordinates to CDS-annotated
#            reads using a strand-aware offset (default 15 nt from
#            the read 3' end to the P-site). Outputs P_cds, A_cds,
#            E_cds positions and in-CDS flags for downstream pause-
#            score analysis.
# Input:     04_monosome/genome-aligned-bam/<sample>_cds_annotation.csv
# Output:    04_monosome/genome-aligned-bam/<sample>_cds_annotation_with_PEA_offset15_fixedStrand.csv
# Parameters: offset = 15 nt (3' end -> P-site), step = 3 nt (codon)
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
missing <- file_paths[!file.exists(file_paths)]
if (length(missing) > 0) {
  stop("Missing file(s):\n", paste(missing, collapse = "\n"),
       "\n\nFiles currently available in ", input_dir, ":\n",
       paste(list.files(input_dir), collapse = "\n"))
}

# =====================================================================
# 1) Parameters
# =====================================================================
offset <- 15L  # 3' end to P-site offset (nt)
step   <- 3L   # Distance between A/P/E sites: one codon (= 3 nt)

# =====================================================================
# 2) Core function: strand-aware P/A/E-site assignment
# =====================================================================
add_PEA_fixed_strand <- function(dt, offset = 15L, step = 3L) {

  # ---- Validate required columns ----
  need_cols <- c("cds_id", "strand", "cds_5", "cds_3", "read_3", "cds_length")
  miss_cols <- setdiff(need_cols, names(dt))
  if (length(miss_cols) > 0) {
    stop("Missing required column(s): ", paste(miss_cols, collapse = ", "))
  }

  # ---- Standardize data types ----
  dt[, `:=`(
    read_3     = as.integer(read_3),
    cds_5      = as.integer(cds_5),
    cds_3      = as.integer(cds_3),
    cds_length = as.integer(cds_length)
  )]

  # ---- Define genomic left and right CDS boundaries ----
  # Use coordinate values rather than naming convention to determine order
  dt[, cds_left  := pmin(cds_5, cds_3)]
  dt[, cds_right := pmax(cds_5, cds_3)]

  # ---- (A) Calculate genomic P/A/E coordinates strand-aware ----
  # Plus strand:  P = read_3 - offset (translation 5'→3' on genomic +)
  # Minus strand: P = read_3 + offset (translation 5'→3' on genomic -)
  dt[, P_genome := fifelse(strand == "+", read_3 - offset, read_3 + offset)]
  dt[, A_genome := fifelse(strand == "+", P_genome + step, P_genome - step)]
  dt[, E_genome := fifelse(strand == "+", P_genome - step, P_genome + step)]

  # ---- (B) Convert genomic coordinates to CDS coordinates ----
  # Plus strand:  CDS coord increases from cds_left
  # Minus strand: CDS coord decreases from cds_right (translation opposite to genome)
  dt[, P_cds := fifelse(
    strand == "+",
    P_genome - cds_left + 1L,
    cds_right - P_genome + 1L
  )]
  dt[, A_cds := fifelse(
    strand == "+",
    A_genome - cds_left + 1L,
    cds_right - A_genome + 1L
  )]
  dt[, E_cds := fifelse(
    strand == "+",
    E_genome - cds_left + 1L,
    cds_right - E_genome + 1L
  )]

  # ---- (C) Flag whether P/A/E sites fall within the CDS ----
  dt[, P_in_cds := (P_cds >= 1L & P_cds <= cds_length)]
  dt[, A_in_cds := (A_cds >= 1L & A_cds <= cds_length)]
  dt[, E_in_cds := (E_cds >= 1L & E_cds <= cds_length)]

  dt[]
}

# =====================================================================
# 3) Batch processing and output
# =====================================================================
for (fp in file_paths) {
  f <- basename(fp)
  message("Reading: ", fp)
  dt <- fread(fp)

  dt2 <- add_PEA_fixed_strand(dt, offset = offset, step = step)

  out_file <- file.path(
    output_dir,
    str_replace(f, "\\.csv$", sprintf("_with_PEA_offset%d_fixedStrand.csv", offset))
  )
  fwrite(dt2, out_file)
  message("Saved: ", out_file)

  # ---- Quick QC: fraction of P/A/E sites within CDS by strand ----
  qc <- dt2[, .(
    n = .N,
    frac_P_in = mean(P_in_cds, na.rm = TRUE),
    frac_A_in = mean(A_in_cds, na.rm = TRUE),
    frac_E_in = mean(E_in_cds, na.rm = TRUE)
  ), by = strand]

  message("QC (", f, "):")
  print(qc)
}

message("All done.")
