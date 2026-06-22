# ============================================================
# Script:    calculate_offset_table.R
# Pipeline:  Step 05 - Monosome (Ribosome Profiling)
# Purpose:   Calculate reading frame periodicity across a range of
#            P-site offsets (10-18 nt) for each read length (26-30 nt).
#            Identifies the optimal offset per read length by selecting
#            the offset that maximizes the percentage of reads in the
#            dominant frame.
# Input:     05_monosome/genome-aligned-bam/<sample>_cds_annotation.csv
# Output:    05_monosome/offset_optimization/
#              offset_optimization_details.csv
#              offset_optimization_summary.csv
#              best_offsets.csv
# Author:    Guchen-Weng
# Date:      2025-01
# ============================================================

# ---- Libraries ----
suppressPackageStartupMessages({
  library(tidyverse)
})

# ---- Paths ----
# Run from Riboseq_analysis/ root directory
input_dir <- "05_monosome/genome-aligned-bam"
output_dir <- "05_monosome/offset_optimization"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ---- Parameters ----
target_lengths <- c(26, 27, 28, 29, 30)
test_offsets <- 10:18

# ---- Load Input Files ----
files <- list.files(input_dir, pattern = "_cds_annotation.csv$", full.names = TRUE)
# Use WT samples as representative set for offset calibration
files <- files[grep("WT", files)]

if (length(files) == 0) stop("No WT input files found.")

cat("Testing offsets", paste(test_offsets, collapse = ","),
    "for lengths", paste(target_lengths, collapse = ","), "...\n")

# ---- Compute Frame Percentages per Offset ----
results_list <- list()

for (f in files) {
  sample_name <- gsub("_cds_annotation.csv", "", basename(f))
  cat("  Processing Sample:", sample_name, "\n")

  data <- read_csv(f, show_col_types = FALSE) %>%
    filter(read_length %in% target_lengths)

  if (nrow(data) == 0) next

  for (L in target_lengths) {
    sub_data <- data %>% filter(read_length == L)
    if (nrow(sub_data) == 0) next

    for (off in test_offsets) {
      # P-site calculation:
      #   + strand: psite = read_3 - offset
      #   - strand: psite = read_3 + offset
      # Frame determination:
      #   + strand: (psite - cds_5) %% 3
      #   - strand: (cds_5 - psite) %% 3

      # Positive strand
      plus_idx <- which(sub_data$strand == "+")
      frame_plus <- ((sub_data$read_3[plus_idx] - off) - sub_data$cds_5[plus_idx]) %% 3

      # Negative strand
      minus_idx <- which(sub_data$strand == "-")
      frame_minus <- (sub_data$cds_5[minus_idx] - (sub_data$read_3[minus_idx] + off)) %% 3

      frames <- c(frame_plus, frame_minus)
      counts <- table(factor(frames, levels = 0:2))

      res_entry <- tibble(
        Sample = sample_name,
        ReadLength = L,
        Offset = off,
        Frame0 = as.integer(counts["0"]),
        Frame1 = as.integer(counts["1"]),
        Frame2 = as.integer(counts["2"]),
        Total = sum(counts)
      ) %>%
        mutate(
          Frame0 = replace_na(Frame0, 0),
          Frame1 = replace_na(Frame1, 0),
          Frame2 = replace_na(Frame2, 0),
          Total = Frame0 + Frame1 + Frame2
        )

      if (res_entry$Total > 0) {
        res_entry <- res_entry %>%
          mutate(
            Pct0 = Frame0 / Total * 100,
            Pct1 = Frame1 / Total * 100,
            Pct2 = Frame2 / Total * 100,
            MaxPct = max(c(Pct0, Pct1, Pct2))
          )
      } else {
        res_entry <- res_entry %>%
          mutate(Pct0 = 0, Pct1 = 0, Pct2 = 0, MaxPct = 0)
      }

      results_list[[paste(sample_name, L, off)]] <- res_entry
    }
  }
}

# ---- Aggregate and Summarize ----
final_results <- bind_rows(results_list)

summary_results <- final_results %>%
  group_by(ReadLength, Offset) %>%
  summarise(
    MeanPct0 = mean(Pct0),
    MeanPct1 = mean(Pct1),
    MeanPct2 = mean(Pct2),
    MeanMaxPct = mean(MaxPct),
    DominantFrame = as.numeric(names(which.max(c(mean(Pct0), mean(Pct1), mean(Pct2))))) - 1,
    .groups = "drop"
  )

# ---- Identify Best Offsets ----
best_offsets <- summary_results %>%
  group_by(ReadLength) %>%
  filter(MeanMaxPct == max(MeanMaxPct)) %>%
  select(ReadLength, Offset, MeanMaxPct, DominantFrame)

# ---- Save Outputs ----
write_csv(final_results, file.path(output_dir, "offset_optimization_details.csv"))
write_csv(summary_results, file.path(output_dir, "offset_optimization_summary.csv"))
write_csv(best_offsets, file.path(output_dir, "best_offsets.csv"))

cat("Analysis complete.\n")
print(best_offsets)
