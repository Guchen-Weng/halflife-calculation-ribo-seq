# ============================================================
# Script:    process_trim_overlap.R
# Pipeline:  Step 05 - Monosome (Ribosome Profiling)
# Purpose:   Trim 20 nt from each CDS end (40 nt total), aggregate
#            read counts per CDS, compute Reads Per Codon (RPC),
#            filter by RPC >= 0.1, and analyze replicate overlap
#            (Jaccard index) for WT and EV conditions.
# Input:     05_monosome/filter-26-30-off13/<sample>_cds_annotation.csv
# Output:    05_monosome/filter-26-30-off13/trim/
#              <sample>_trimmed.csv
#              overlap_stats.txt
# Author:    Guchen-Weng
# Date:      2025-01
# ============================================================

# ---- Paths ----
# Run from Riboseq_analysis/ root directory
input_dir <- "05_monosome/filter-26-30-off13"
output_dir <- "05_monosome/filter-26-30-off13/trim"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ---- List Input Files ----
files <- list.files(input_dir, pattern = "_cds_annotation.csv$", full.names = TRUE)

# ---- Process a Single File ----
process_file <- function(file_path) {
  message(paste("Processing", basename(file_path), "..."))

  data <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)

  # Extract unique CDS info (cds_length is constant per CDS)
  cds_info <- unique(data[, c("cds_id", "cds_length")])

  # Aggregate read counts per CDS
  read_counts <- aggregate(readcount ~ cds_id, data = data, sum)

  # Merge CDS info with aggregated counts
  merged <- merge(cds_info, read_counts, by = "cds_id")

  # Trim 20 nt from each CDS end (40 nt total)
  merged$cds_length_trimmed <- merged$cds_length - 40

  # Remove CDS too short after trimming
  merged <- merged[merged$cds_length_trimmed > 0, ]

  # Calculate Reads Per Codon (RPC)
  # RPC = sum_read_count * 3 / cds_length_trimmed
  merged$readpercodon <- (merged$readcount * 3) / merged$cds_length_trimmed

  # Filter: RPC >= 0.1
  final_data <- merged[merged$readpercodon >= 0.1, ]

  # Select and rename output columns
  output_data <- final_data[, c("cds_id", "cds_length_trimmed", "readcount", "readpercodon")]
  colnames(output_data) <- c("cds_id", "cds_length_trimmed", "read_count", "readpercodon")

  # Write trimmed CSV
  base_name <- basename(file_path)
  sample_name <- sub("_cds_annotation.csv", "", base_name)
  output_filename <- paste0(sample_name, "_trimmed.csv")
  output_path <- file.path(output_dir, output_filename)
  write.csv(output_data, output_path, row.names = FALSE, quote = FALSE)

  return(list(sample = sample_name, ids = output_data$cds_id))
}

# ---- Process All Files ----
results <- list()
for (f in files) {
  res <- process_file(f)
  results[[res$sample]] <- res$ids
}

# ---- Overlap Analysis ----
message("\nPerforming Overlap Analysis...")

calc_overlap <- function(set1, set2, name1, name2) {
  common <- intersect(set1, set2)
  u1 <- setdiff(set1, set2)
  u2 <- setdiff(set2, set1)

  n_common <- length(common)
  n_total <- length(union(set1, set2))
  pct_overlap <- if (n_total > 0) (n_common / n_total) * 100 else 0

  cat(sprintf("\nOverlap: %s vs %s\n", name1, name2))
  cat(sprintf("  %s unique: %d\n", name1, length(u1)))
  cat(sprintf("  %s unique: %d\n", name2, length(u2)))
  cat(sprintf("  Common: %d\n", n_common))
  cat(sprintf("  Jaccard Index (Overlap Rate): %.2f%%\n", pct_overlap))

  return(common)
}

# Within-group overlap
wt_overlap <- calc_overlap(results[["WT_1"]], results[["WT_2"]], "WT_1", "WT_2")
ev_overlap <- calc_overlap(results[["EV_1"]], results[["EV_2"]], "EV_1", "EV_2")

# Cross-group overlap of consensus CDS sets
final_overlap <- calc_overlap(wt_overlap, ev_overlap, "WT_Overlap", "EV_Overlap")

# ---- Save Overlap Statistics ----
sink(file.path(output_dir, "overlap_stats.txt"))
cat("Overlap Analysis Results\n")
cat("========================\n")

cat(sprintf("\n--- WT_1 vs WT_2 ---\n"))
common_wt <- intersect(results[["WT_1"]], results[["WT_2"]])
cat(sprintf("WT_1 Count: %d\n", length(results[["WT_1"]])))
cat(sprintf("WT_2 Count: %d\n", length(results[["WT_2"]])))
cat(sprintf("Intersection: %d\n", length(common_wt)))
cat(sprintf("Union: %d\n", length(union(results[["WT_1"]], results[["WT_2"]]))))
cat(sprintf("Overlap Rate (Jaccard): %.2f%%\n",
    length(common_wt) / length(union(results[["WT_1"]], results[["WT_2"]])) * 100))

cat(sprintf("\n--- EV_1 vs EV_2 ---\n"))
common_ev <- intersect(results[["EV_1"]], results[["EV_2"]])
cat(sprintf("EV_1 Count: %d\n", length(results[["EV_1"]])))
cat(sprintf("EV_2 Count: %d\n", length(results[["EV_2"]])))
cat(sprintf("Intersection: %d\n", length(common_ev)))
cat(sprintf("Union: %d\n", length(union(results[["EV_1"]], results[["EV_2"]]))))
cat(sprintf("Overlap Rate (Jaccard): %.2f%%\n",
    length(common_ev) / length(union(results[["EV_1"]], results[["EV_2"]])) * 100))

cat(sprintf("\n--- WT_Overlap vs EV_Overlap ---\n"))
final_common <- intersect(common_wt, common_ev)
cat(sprintf("WT_Overlap Count: %d\n", length(common_wt)))
cat(sprintf("EV_Overlap Count: %d\n", length(common_ev)))
cat(sprintf("Intersection: %d\n", length(final_common)))
cat(sprintf("Union: %d\n", length(union(common_wt, common_ev))))
cat(sprintf("Overlap Rate (Jaccard): %.2f%%\n",
    length(final_common) / length(union(common_wt, common_ev)) * 100))

sink()

message(paste("\nAnalysis complete. Results saved to", output_dir))
