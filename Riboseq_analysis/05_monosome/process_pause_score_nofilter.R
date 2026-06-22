# ============================================================
# Script:    process_pause_score_nofilter.R
# Pipeline:  Step 05 - Monosome (Ribosome Profiling)
# Purpose:   Compute ribosome pause scores from filtered P-site
#            annotations. For each observed P-site position, calculate
#            the local pause score (observed counts / gene mean density).
#            Extract codons at P-sites and summarize mean pause scores
#            per codon. Perform WT vs EV comparison with log2 fold change.
# Input:     05_monosome/filter-26-30-off13/<sample>_cds_annotation.csv
#            05_monosome/filter-26-30-off13/trim/<sample>_trimmed.csv
#            reference/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.dna.chromosome.Chromosome.fa
# Output:    05_monosome/filter-26-30-off13/pausescore-0nofilter/
#              <sample>_filtered_annotation.csv
#              <sample>_filter_stats.txt
#              <sample>_codon_summary.csv
#              WT_vs_EV_pause_score_comparison.csv
# Author:    Guchen-Weng
# Date:      2025-01
# ============================================================

# ---- Libraries ----
suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# ---- Paths ----
# Run from Riboseq_analysis/ root directory
input_dir <- "05_monosome/filter-26-30-off13"
trim_dir <- "05_monosome/filter-26-30-off13/trim"
output_dir <- "05_monosome/filter-26-30-off13/pausescore-0nofilter"
fasta_file <- "reference/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.dna.chromosome.Chromosome.fa"
samples <- c("WT_1", "WT_2", "EV_1", "EV_2")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ---- Load Genome Sequence ----
cat("Loading genome sequence...\n")
if (!file.exists(fasta_file)) stop(paste("FASTA file not found:", fasta_file))
genome <- readDNAStringSet(fasta_file)
genome_seq <- genome[[1]]

# ---- Helper: Vectorized Codon Extraction ----
get_codons_vectorized <- function(psites, strands, genome_seq) {
  n <- length(psites)
  codons <- character(n)

  # Positive strand: read forward from psite
  pos_idx <- which(strands == "+")
  if (length(pos_idx) > 0) {
    starts <- psites[pos_idx]
    ends <- starts + 2
    valid <- ends <= length(genome_seq)

    if (any(valid)) {
      v <- Views(genome_seq, start = starts[valid], end = ends[valid])
      codons[pos_idx[valid]] <- as.character(v)
    }
    if (any(!valid)) codons[pos_idx[!valid]] <- NA
  }

  # Negative strand: read reverse complement
  neg_idx <- which(strands == "-")
  if (length(neg_idx) > 0) {
    ends <- psites[neg_idx]
    starts <- ends - 2
    valid <- starts >= 1

    if (any(valid)) {
      v <- Views(genome_seq, start = starts[valid], end = ends[valid])
      codons[neg_idx[valid]] <- as.character(reverseComplement(v))
    }
    if (any(!valid)) codons[neg_idx[!valid]] <- NA
  }

  return(codons)
}

# ---- Process a Single Sample ----
process_sample <- function(sample_name) {
  cat(sprintf("Processing %s...\n", sample_name))

  # 1. Load files
  annot_file <- file.path(input_dir, paste0(sample_name, "_cds_annotation.csv"))
  trimmed_file <- file.path(trim_dir, paste0(sample_name, "_trimmed.csv"))

  if (!file.exists(annot_file) || !file.exists(trimmed_file)) {
    warning(paste("Files not found for", sample_name))
    return(NULL)
  }

  annot_data <- read.csv(annot_file, stringsAsFactors = FALSE)
  trimmed_data <- read.csv(trimmed_file, stringsAsFactors = FALSE)

  # 2. Filter annotation to CDS that passed RPC threshold
  valid_ids <- unique(trimmed_data$cds_id)

  n_total <- nrow(annot_data)
  filtered_annot <- annot_data %>% filter(cds_id %in% valid_ids)
  n_kept <- nrow(filtered_annot)
  n_removed <- n_total - n_kept

  cat(sprintf("  Filtered CDS IDs: Kept %d / %d rows (Removed %d)\n",
              n_kept, n_total, n_removed))

  # Save filter statistics
  stats_file <- file.path(output_dir, paste0(sample_name, "_filter_stats.txt"))
  writeLines(c(
    paste("Original rows:", n_total),
    paste("Filtered rows:", n_kept),
    paste("Removed rows:", n_removed),
    paste("Unique CDS in Trimmed:", length(valid_ids)),
    paste("Unique CDS in Filtered Annot:", length(unique(filtered_annot$cds_id)))
  ), stats_file)

  # Save filtered annotation
  filtered_csv <- file.path(output_dir, paste0(sample_name, "_filtered_annotation.csv"))
  write.csv(filtered_annot, filtered_csv, row.names = FALSE)

  # 3. Calculate Pause Scores (no zero-filling: only observed positions)
  cat("  Calculating Pause Scores (no zero filling)...\n")

  # 3a: Aggregate reads to P-site positions
  position_counts <- filtered_annot %>%
    group_by(cds_id, psite, strand) %>%
    summarise(C_ij = sum(readcount), .groups = "drop")

  # 3b: Calculate per-gene mean density
  gene_meta <- trimmed_data %>%
    select(cds_id, cds_length_trimmed) %>%
    distinct()

  gene_stats <- position_counts %>%
    group_by(cds_id) %>%
    summarise(
      internal_total_counts = sum(C_ij),
      .groups = "drop"
    ) %>%
    inner_join(gene_meta, by = "cds_id") %>%
    mutate(
      internal_nt = cds_length_trimmed,
      mean_density_nt = internal_total_counts / internal_nt
    )

  # 3c: Compute local pause scores
  # PS_ij = C_ij / mean_density_nt (for observed positions only)
  data_with_score <- position_counts %>%
    inner_join(gene_stats, by = "cds_id") %>%
    mutate(
      pause_score = ifelse(mean_density_nt > 0, C_ij / mean_density_nt, 0)
    )

  # 4. Extract codons at P-sites
  cat("  Extracting codons for observed positions...\n")
  data_with_score$codon <- get_codons_vectorized(
    data_with_score$psite, data_with_score$strand, genome_seq
  )

  # Remove invalid codons (NA, wrong length, ambiguous bases)
  data_with_score <- data_with_score %>%
    filter(!is.na(codon), nchar(codon) == 3, !grepl("N", codon))

  # 5. Summarize pause scores by codon
  cat("  Summarizing by codon...\n")
  codon_summary <- data_with_score %>%
    group_by(codon) %>%
    summarise(
      Mean_pause = mean(pause_score),
      N_instances = n(),
      Total_counts = sum(C_ij)
    ) %>%
    mutate(Sample = sample_name)

  # Save codon summary
  write.csv(codon_summary,
            file.path(output_dir, paste0(sample_name, "_codon_summary.csv")),
            row.names = FALSE)

  return(codon_summary)
}

# ---- Run Analysis for All Samples ----
results_list <- list()

for (sample in samples) {
  res <- process_sample(sample)
  if (!is.null(res)) {
    results_list[[sample]] <- res
  }
}

# ---- WT vs EV Comparison ----
cat("Performing WT vs EV comparison...\n")

if (length(results_list) == 4) {
  all_results <- bind_rows(results_list)

  # Add group label
  all_results <- all_results %>%
    mutate(Group = ifelse(grepl("WT", Sample), "WT", "EV"))

  # Average pause score per codon per group
  grouped_stats <- all_results %>%
    group_by(Group, codon) %>%
    summarise(
      Avg_Mean_Pause = mean(Mean_pause),
      .groups = "drop"
    )

  # Pivot to wide format: WT and EV columns
  wide_stats <- grouped_stats %>%
    tidyr::pivot_wider(names_from = Group, values_from = Avg_Mean_Pause)

  # Calculate log2 fold change (WT / EV) with small pseudocount
  epsilon <- 1e-5
  wide_stats <- wide_stats %>%
    mutate(
      log2FC = log2((WT + epsilon) / (EV + epsilon))
    ) %>%
    arrange(desc(log2FC))

  # Save comparison table
  write.csv(wide_stats,
            file.path(output_dir, "WT_vs_EV_pause_score_comparison.csv"),
            row.names = FALSE)
  cat("Comparison saved to WT_vs_EV_pause_score_comparison.csv\n")

} else {
  warning("Not all samples processed successfully. Skipping comparison.")
}

cat("Done.\n")
