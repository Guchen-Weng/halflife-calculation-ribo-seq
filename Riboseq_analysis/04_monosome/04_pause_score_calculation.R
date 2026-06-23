#!/usr/bin/env Rscript
# ============================================================
# Script:    04_pause_score_calculation.R
# Pipeline:  Step 04 - Monosome (Ribosome Profiling)
# Purpose:   Calculate codon-level pause scores from CDS-annotated
#            read tables with pre-computed P/A/E-site positions
#            (offset 15). Computes per-codon, per-amino-acid, and
#            gene-half asymmetry scores for each of 4 samples
#            (WT_1, WT_2, EV_1, EV_2).
# Input:     *_cds_annotation_with_PEA_offset15_fixedStrand.csv
#            ../reference/Escherichia_coli_...cds.all.fa
# Output:    PauseScore_fixedStrand_len26-32/
#              *_site_codon_nt_table.csv
#              *_pause_scores_codon_A_P_E.csv
#              *_pause_scores_aa_A_P_E.csv
#              *_asymmetry_scores_gene.csv
# Parameters: read length 26-32, exclude first 27nt / last 12nt,
#             gene-level min density 0.5 reads/nt
# Author:    Guchen-Weng
# Date:      2025-01
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(Biostrings)
  library(stringr)
})

# =====================================================================
# 0) Paths and input files
# =====================================================================
# Run from Riboseq_analysis/ root directory
input_dir  <- "04_monosome/genome-aligned-bam"
output_dir <- "04_monosome/PauseScore_fixedStrand_len26-32"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

files <- file.path(input_dir, c(
  "WT_1_cds_annotation_with_PEA_offset15_fixedStrand.csv",
  "WT_2_cds_annotation_with_PEA_offset15_fixedStrand.csv",
  "EV_1_cds_annotation_with_PEA_offset15_fixedStrand.csv",
  "EV_2_cds_annotation_with_PEA_offset15_fixedStrand.csv"
))
stopifnot(all(file.exists(files)))

cds_fa <- "reference/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.cds.all.fa"
stopifnot(file.exists(cds_fa))

# =====================================================================
# 1) Parameters
# =====================================================================
len_min <- 26L
len_max <- 32L

# Exclude terminal CDS regions
exclude_after_start_nt <- 27L
exclude_before_stop_nt <- 12L

# Gene-level filter: discard genes with mean density < 0.5 reads/nt
min_reads_per_nt <- 0.5

window_len <- function(L) {
  (L - exclude_before_stop_nt) - (1L + exclude_after_start_nt) + 1L
}

# =====================================================================
# 2) Read CDS FASTA and build codon map
# =====================================================================
cds <- readDNAStringSet(cds_fa)
cds_dt <- data.table(
  cds_id_raw = names(cds),
  seq = as.character(cds)
)

# Standardize IDs to CDS:xxxx format
cds_dt[, cds_id := sub(" .*", "", cds_id_raw)]
cds_dt[, cds_id := ifelse(str_detect(cds_id, "^CDS:"),
                          cds_id, paste0("CDS:", cds_id))]
cds_dt[, cds_length_fa := nchar(seq)]
setkey(cds_dt, cds_id)

make_codon_map <- function(cds_id, seq) {
  L <- nchar(seq)
  ncod <- L %/% 3
  starts <- seq(1, by = 3, length.out = ncod)
  codons <- substring(seq, starts, starts + 2)
  data.table(cds_id = cds_id, codon_idx = seq_len(ncod), codon_nt = codons)
}
codon_map <- rbindlist(Map(make_codon_map, cds_dt$cds_id, cds_dt$seq))
setkey(codon_map, cds_id, codon_idx)

aa_table <- Biostrings::GENETIC_CODE

# =====================================================================
# 3) Process one sample: A/P/E-site nucleotide-density profiles
#    Uses pre-computed P_cds/A_cds/E_cds coordinates from input table
# =====================================================================
calc_pause_one_sample <- function(f) {
  dt <- fread(f)

  need_cols <- c("cds_id", "strand", "read_length", "readcount",
                 "P_cds", "A_cds", "E_cds",
                 "P_in_cds", "A_in_cds", "E_in_cds")
  miss <- setdiff(need_cols, names(dt))
  if (length(miss) > 0) {
    stop("Missing required columns in: ", f,
         "\nMissing: ", paste(miss, collapse = ", "))
  }

  dt[, `:=`(
    read_length = as.integer(read_length),
    readcount   = as.numeric(readcount)
  )]

  # (A) Retain reads of length 26-32 nt
  dt <- dt[read_length >= len_min & read_length <= len_max]
  if (nrow(dt) == 0) stop("No data after read-length filter: ", f)

  # (B) Add CDS lengths from FASTA
  dt <- cds_dt[dt, on = "cds_id", nomatch = 0]
  if (nrow(dt) == 0) stop("No cds_id matched to FASTA: ", f)

  # (C) Keep only in-CDS sites for pause-score analysis
  make_site_dt <- function(site_name, pos_col, in_col) {
    x <- dt[get(in_col) == TRUE]
    x[, site := site_name]
    x[, site_pos_cds := as.integer(get(pos_col))]
    x[, .(site, cds_id, cds_length_fa, site_pos_cds, readcount)]
  }

  dA <- make_site_dt("A", "A_cds", "A_in_cds")
  dP <- make_site_dt("P", "P_cds", "P_in_cds")
  dE <- make_site_dt("E", "E_cds", "E_in_cds")
  dtx <- rbindlist(list(dA, dP, dE), use.names = TRUE)

  # (D) Exclude CDS termini
  dtx <- dtx[
    site_pos_cds >= (1L + exclude_after_start_nt) &
    site_pos_cds <= (cds_length_fa - exclude_before_stop_nt)
  ]
  dtx <- dtx[window_len(cds_length_fa) > 0]
  if (nrow(dtx) == 0) stop("No data after terminal filter: ", f)

  # (E) Aggregate per-nucleotide coverage
  nt_den <- dtx[, .(cov = sum(readcount, na.rm = TRUE)),
                by = .(site, cds_id, cds_length_fa, site_pos_cds)]

  # (F) Gene-level filter: >= 0.5 reads/nt
  gene_cov <- nt_den[, .(total = sum(cov)),
                     by = .(site, cds_id, cds_length_fa)]
  gene_cov[, win_nt := window_len(cds_length_fa)]
  gene_cov[, reads_per_nt := total / win_nt]

  keep <- gene_cov[reads_per_nt >= min_reads_per_nt, .(site, cds_id)]
  setkey(keep, site, cds_id)
  setkey(nt_den, site, cds_id)
  nt_den <- keep[nt_den, nomatch = 0]
  if (nrow(nt_den) == 0) stop("No data after gene filter: ", f)

  # (G) Pause-score normalization: coverage / mean coverage per gene
  nt_den[, win_nt := window_len(cds_length_fa)]
  nt_den[, mean_cov := sum(cov) / win_nt[1], by = .(site, cds_id)]
  nt_den[, norm := fifelse(mean_cov > 0, cov / mean_cov, NA_real_)]

  # (H) Map nucleotides to codon indices
  nt_den[, codon_idx := as.integer((site_pos_cds + 2L) %/% 3L)]
  codon_den <- nt_den[, .(
    pause = mean(norm, na.rm = TRUE),
    n_nt = .N
  ), by = .(site, cds_id, codon_idx)]

  # (I) Add codon sequence and amino acid annotation
  codon_den <- codon_map[codon_den, on = .(cds_id, codon_idx)]
  codon_den <- codon_den[!is.na(codon_nt)]

  codon_den[, `:=`(
    aa = aa_table[codon_nt],
    codon_nt1 = substr(codon_nt, 1, 1),
    codon_nt2 = substr(codon_nt, 2, 2),
    codon_nt3 = substr(codon_nt, 3, 3)
  )]

  # Table 1: Per-gene, per-position
  site_codon_nt_table <- codon_den[, .(
    site, cds_id, codon_idx, codon_nt,
    codon_nt1, codon_nt2, codon_nt3,
    aa, pause, n_nt
  )]

  # Table 2: Mean pause score per codon (64 codons, A/P/E sites)
  codon_pause <- codon_den[, .(
    pause_mean = mean(pause, na.rm = TRUE),
    n = .N
  ), by = .(site, codon = codon_nt)][order(site, -pause_mean)]
  codon_pause[, aa := aa_table[codon]]

  # Table 3: Amino-acid-level mean pause scores
  aa_pause <- codon_pause[!is.na(aa),
    .(pause_mean = weighted.mean(pause_mean, w = n),
      n_codons = sum(n)),
    by = .(site, aa)][order(site, -pause_mean)]

  # Table 4: P-site coverage asymmetry (first vs. second half of gene)
  P_nt <- nt_den[site == "P"]
  P_nt[, half := fifelse(site_pos_cds <= (cds_length_fa %/% 2L),
                          "first", "second")]
  asym <- P_nt[, .(
    d1 = sum(cov[half == "first"], na.rm = TRUE),
    d2 = sum(cov[half == "second"], na.rm = TRUE)
  ), by = .(cds_id)]
  asym[, asymmetry_log2 := log2((d2 + 1e-9) / (d1 + 1e-9))]

  sample_name <- str_replace(basename(f), "\\.csv$", "")
  list(sample = sample_name,
       site_codon_nt_table = site_codon_nt_table,
       codon_pause = codon_pause,
       aa_pause = aa_pause,
       asym = asym)
}

# =====================================================================
# 4) Process all four samples and write output
# =====================================================================
for (f in files) {
  res <- calc_pause_one_sample(f)

  fwrite(res$site_codon_nt_table,
         file.path(output_dir, paste0(res$sample, "_site_codon_nt_table.csv")))
  fwrite(res$codon_pause,
         file.path(output_dir, paste0(res$sample, "_pause_scores_codon_A_P_E.csv")))
  fwrite(res$aa_pause,
         file.path(output_dir, paste0(res$sample, "_pause_scores_aa_A_P_E.csv")))
  fwrite(res$asym,
         file.path(output_dir, paste0(res$sample, "_asymmetry_scores_gene.csv")))

  message("Finished: ", res$sample)
}

message("All done. Output in: ", output_dir)
