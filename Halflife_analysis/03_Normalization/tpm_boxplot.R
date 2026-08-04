# ============================================================
# Script:    tpm_boxplot.R
# Pipeline:  Step 03 - Normalization
# Purpose:   Calculate TPM from raw counts and generate a boxplot
#            of log2(TPM+1) distribution across all samples
# Input:     ../../02-Readcount/merged_gene_counts_final.csv
# Output:    ../output/tpm_boxplot_continuous.pdf
#            ../output/tpm_log_data_long.csv
#            ../output/tpm_log_matrix.csv
# Author:    Guchen-Weng
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)

# Input file is two levels up from 03_Normalization/
input_file <- "../../02-Readcount/merged_gene_counts_final.csv"
output_dir <- "../output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Read count data
data <- read.csv(input_file, header = TRUE, check.names = FALSE)

# Calculate gene lengths
if (!all(c("start", "end") %in% colnames(data))) {
  stop("Input file must contain 'start' and 'end' columns for TPM calculation.")
}
gene_lengths <- data$end - data$start + 1
gene_lengths[gene_lengths < 1] <- 1

# Identify sample (count) columns
non_sample_cols <- c("Geneid", "chromosome", "start", "end", "strand",
                     "biotype", "description", "gene_name")
sample_cols <- setdiff(colnames(data), non_sample_cols)

# Extract count matrix and fix "omin" -> "0min" typo
counts_matrix <- data[, sample_cols]
colnames(counts_matrix) <- gsub("omin", "0min", colnames(counts_matrix))

# TPM calculation
calculate_tpm <- function(counts, lengths) {
  len_kb <- lengths / 1000
  rpk <- counts / len_kb
  scaling_factors <- colSums(rpk) / 1e6
  scaling_factors[scaling_factors == 0] <- 1
  tpm <- t(t(rpk) / scaling_factors)
  return(tpm)
}

tpm_matrix <- calculate_tpm(counts_matrix, gene_lengths)

# Log2 transformation
log_tpm_matrix <- log2(tpm_matrix + 1)

# Convert to long format for ggplot
plot_data_long <- gather(as.data.frame(log_tpm_matrix),
                         key = "Sample", value = "LogTPM")

# Parse sample names (format: Replicate_Condition_Time_ID)
parse_sample_info <- function(sample_name) {
  parts <- strsplit(sample_name, "_")[[1]]
  condition <- parts[2]
  time_point <- parts[3]
  return(c(Condition = condition, Time = time_point))
}

sample_info <- t(sapply(plot_data_long$Sample, parse_sample_info))
plot_data_long$Condition <- sample_info[, "Condition"]
plot_data_long$Time <- sample_info[, "Time"]

# Rename 2min -> 3min (data correction; use direct replacement to
# avoid accidentally changing 12min -> 13min)
plot_data_long$Time[plot_data_long$Time == "2min"] <- "3min"

# Set factor levels
time_levels <- c("0min", "3min", "6min", "12min")
condition_levels <- c("WT", "EV")

plot_data_long$Time <- factor(plot_data_long$Time, levels = time_levels)
plot_data_long$Condition <- factor(plot_data_long$Condition,
                                   levels = condition_levels)

# Order samples: Time -> Condition -> Sample
unique_samples <- unique(plot_data_long[, c("Sample", "Time", "Condition")])
unique_samples <- unique_samples[order(unique_samples$Time,
                                       unique_samples$Condition,
                                       unique_samples$Sample), ]

# Assign continuous x-axis positions
unique_samples$X_Pos <- 1:nrow(unique_samples)

# Merge positions back
plot_data_long <- merge(plot_data_long,
                        unique_samples[, c("Sample", "X_Pos")],
                        by = "Sample")

# Save intermediate data
write.csv(plot_data_long,
          file = file.path(output_dir, "tpm_log_data_long.csv"),
          row.names = FALSE)
write.csv(as.data.frame(log_tpm_matrix),
          file = file.path(output_dir, "tpm_log_matrix.csv"),
          row.names = TRUE)

# --- Plot ---
colors <- c("WT" = "#d74596", "EV" = "#95c94a")

# Compute x-axis tick positions (center of each time group)
time_breaks <- unique_samples %>%
  group_by(Time) %>%
  summarise(Center = mean(X_Pos))

print("Time Breaks:")
print(time_breaks)

p <- ggplot(plot_data_long,
            aes(x = X_Pos, y = LogTPM, group = X_Pos, fill = Condition)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  scale_fill_manual(values = colors) +
  scale_x_continuous(breaks = time_breaks$Center,
                     labels = as.character(time_breaks$Time)) +
  scale_y_continuous(limits = c(-3, 20), breaks = seq(0, 20, by = 2)) +
  labs(x = "Time", y = "Log2(TPM + 1)",
       title = "TPM Distribution per Sample") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.position = "top"
  )

# Save plot (aspect ratio 1:1.5)
output_plot <- file.path(output_dir, "tpm_boxplot_continuous.pdf")
ggsave(output_plot, plot = p, width = 6, height = 9)
print(paste("Saved plot to:", normalizePath(output_plot)))
