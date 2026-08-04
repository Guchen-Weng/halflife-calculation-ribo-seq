# ============================================================
# Script:    sample_correlation_heatmap.R
# Pipeline:  DESeq2 - Sample Quality Assessment
# Purpose:   Compute sample-to-sample Pearson correlation on
#            log2(TPM+1) values and plot a heatmap with Group
#            (WT/EV) and Time annotations (pheatmap).
# Input:     ../../02-Readcount/merged_gene_counts_final.csv
# Output:    sample_correlation_heatmap.pdf / .png
# Author:    Guchen-Weng
# ============================================================

suppressPackageStartupMessages({
  library(pheatmap)
  library(RColorBrewer)
  library(ggplot2)
})

# ---- Paths ----
input_file <- "../../02-Readcount/merged_gene_counts_final.csv"
output_dir <- "."

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

setwd(output_dir)

# ---- 1. Read data ----
message("Reading input file...")
data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)

# ---- 2. Calculate TPM ----
if (!all(c("start", "end") %in% colnames(data))) {
  stop("Input file must contain 'start' and 'end' columns for TPM calculation.")
}

gene_lengths <- data$end - data$start + 1
gene_lengths[gene_lengths < 1] <- 1

# Identify sample columns
non_sample_cols <- c("Geneid", "chromosome", "start", "end", "strand",
                     "biotype", "description", "gene_name")
sample_cols <- setdiff(colnames(data), non_sample_cols)

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

message("Calculating TPM...")
tpm_matrix <- calculate_tpm(counts_matrix, gene_lengths)

# ---- 3. Aggregate TPM by gene_name ----
message("Aggregating TPM by gene_name...")
gene_names <- data$gene_name
tpm_df <- as.data.frame(tpm_matrix)
tpm_df$gene_name <- gene_names

tpm_agg <- aggregate(. ~ gene_name, data = tpm_df, FUN = sum)
rownames(tpm_agg) <- tpm_agg$gene_name
tpm_agg$gene_name <- NULL
tpm_final <- as.matrix(tpm_agg)

# Remove genes with zero total TPM
tpm_final <- tpm_final[rowSums(tpm_final) > 0, ]

# ---- 4. Parse sample metadata and sort ----
sample_names <- colnames(tpm_final)
meta <- data.frame(sample_name = sample_names, stringsAsFactors = FALSE)

# Expected format: SampleID_Group_Time_Rep (e.g., A1_WT_0min_1)
parts <- strsplit(sample_names, "_")
meta$Group <- sapply(parts, function(x) x[2])
meta$Time <- sapply(parts, function(x) x[3])
meta$Rep <- sapply(parts, function(x) x[4])
meta$TimeNum <- as.numeric(gsub("min", "", meta$Time))

meta$Group <- factor(meta$Group, levels = c("WT", "EV"))

# Order: WT first, then EV; within group by Time
meta_ordered <- meta[order(meta$Group, meta$TimeNum, meta$sample_name), ]

# Reorder TPM matrix
tpm_sorted <- tpm_final[, meta_ordered$sample_name]

# ---- 5. Compute Pearson correlation on log2(TPM+1) ----
message("Calculating correlation matrix...")
log_tpm <- log2(tpm_sorted + 1)
cor_mat <- cor(log_tpm, method = "pearson")

# ---- 6. Plot heatmap with annotations ----
annotation_col <- data.frame(
  Group = meta_ordered$Group,
  Time = factor(meta_ordered$Time,
                levels = unique(meta_ordered$Time[order(meta_ordered$TimeNum)]))
)
rownames(annotation_col) <- meta_ordered$sample_name

ann_colors <- list(
  Group = c(WT = "#1B9E77", EV = "#D95F02"),
  Time = c("0min"  = "#E41A1C", "2min"  = "#377EB8",
           "6min"  = "#4DAF4A", "12min" = "#984EA3")
)

message("Generating heatmap...")
pheatmap(cor_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_colnames = TRUE,
         show_rownames = TRUE,
         main = "Sample Correlation Heatmap (Pearson, log2(TPM+1))",
         filename = "sample_correlation_heatmap.pdf",
         width = 10,
         height = 8)

pheatmap(cor_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_colnames = TRUE,
         show_rownames = TRUE,
         main = "Sample Correlation Heatmap (Pearson, log2(TPM+1))",
         filename = "sample_correlation_heatmap.png",
         width = 10,
         height = 8)

message("Done! Heatmaps saved in ", output_dir)
