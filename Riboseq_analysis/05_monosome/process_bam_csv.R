# ============================================================
# Script:    process_bam_csv.R
# Pipeline:  Step 05 - Monosome (Ribosome Profiling)
# Purpose:   Filter annotated reads by read length (26-30 nt),
#            calculate P-site position (offset 13), and remove reads
#            whose P-site falls within 20 nt of either CDS end.
# Input:     05_monosome/genome-aligned-bam/<sample>_cds_annotation.csv
# Output:    05_monosome/filter-26-30-off13/<sample>_cds_annotation.csv
# Author:    Guchen-Weng
# Date:      2025-01
# ============================================================

# ---- Paths ----
# Run from Riboseq_analysis/ root directory
input_dir <- "05_monosome/genome-aligned-bam"
output_dir <- "05_monosome/filter-26-30-off13"

# Create output directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat(sprintf("Created directory: %s\n", output_dir))
}

# ---- List Input CSV Files ----
files <- list.files(path = input_dir, pattern = "\\.csv$", full.names = TRUE)

if (length(files) == 0) {
  stop("No CSV files found in ", input_dir, " directory.")
}

cat(sprintf("Found %d CSV files to process.\n", length(files)))

# ---- Process Each File ----
for (file_path in files) {
  filename <- basename(file_path)
  cat(sprintf("------------------------------\n"))
  cat(sprintf("Processing %s...\n", filename))

  data <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  original_count <- nrow(data)

  # ---- Filter 1: Read length 26-30 nt ----
  data_filtered <- data[data$read_length >= 26 & data$read_length <= 30, ]
  count_after_len <- nrow(data_filtered)

  if (count_after_len == 0) {
    cat(sprintf("  Processed: %d\n  Kept: 0\n  Filtered (length): %d\n  Filtered (position): 0\n",
                original_count, original_count))
    next
  }

  # ---- Calculate P-site (Offset 13) ----
  # + strand: psite = read_3 - 13
  # - strand: psite = read_3 + 13
  psite <- numeric(nrow(data_filtered))

  plus_indices <- which(data_filtered$strand == "+")
  minus_indices <- which(data_filtered$strand == "-")

  if (length(plus_indices) > 0) {
    psite[plus_indices] <- data_filtered$read_3[plus_indices] - 13
  }
  if (length(minus_indices) > 0) {
    psite[minus_indices] <- data_filtered$read_3[minus_indices] + 13
  }

  data_filtered$psite <- psite

  # ---- Filter 2: P-site distance to CDS ends > 20 nt ----
  keep_indices <- which(abs(data_filtered$psite - data_filtered$cds_5) > 20 &
                        abs(data_filtered$psite - data_filtered$cds_3) > 20)

  final_data <- data_filtered[keep_indices, ]
  final_count <- nrow(final_data)

  # ---- Save Output ----
  output_path <- file.path(output_dir, filename)
  write.csv(final_data, output_path, row.names = FALSE, quote = FALSE)

  cat(sprintf("  Processed: %d\n", original_count))
  cat(sprintf("  Kept: %d\n", final_count))
  cat(sprintf("  Filtered (length): %d\n", original_count - count_after_len))
  cat(sprintf("  Filtered (position): %d\n", count_after_len - final_count))
}

cat("------------------------------\n")
cat("All files processed successfully.\n")
