# ============================================================
# Script:    run_deseq2.R
# Pipeline:  DESeq2 - Differential Expression Analysis
# Purpose:   Run DESeq2 comparing WT vs EV at the 0min time point.
#            Aggregates counts by gene_name, generates volcano plot
#            (ggrepel-labeled top 10 up/down genes), and saves full
#            results table.
# Input:     ../../02-Readcount/merged_gene_counts_final.csv
# Output:    deseq2_results_0min.csv
#            volcano_plot_0min.pdf / .png
# Author:    Guchen-Weng
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
})

# ---- Paths ----
input_file <- "../../02-Readcount/merged_gene_counts_final.csv"
output_dir <- "."

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

setwd(output_dir)

# ---- 1. Read and prepare data ----
message("Reading input file...")
data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)

# Use gene_name as identifier; count columns start from column 9
message("Processing gene names...")
counts_raw <- data[, 9:ncol(data)]
gene_names <- data$gene_name

# Handle duplicate gene names by summing counts
if (any(duplicated(gene_names))) {
  message("Duplicate gene names found. Aggregating counts by sum...")
  counts_raw$gene_name <- gene_names
  counts_agg <- aggregate(. ~ gene_name, data = counts_raw, FUN = sum)
  rownames(counts_agg) <- counts_agg$gene_name
  counts_agg$gene_name <- NULL
  counts <- counts_agg
} else {
  rownames(counts_raw) <- gene_names
  counts <- counts_raw
}

# Fix column name typo: "omin" -> "0min"
colnames(counts) <- gsub("omin", "0min", colnames(counts))

# ---- 2. Build sample metadata ----
sample_names <- colnames(counts)
meta <- data.frame(sample_name = sample_names)

# Expected format: SampleID_Group_Time_Rep (e.g., A1_WT_0min_1)
parts <- strsplit(sample_names, "_")
meta$condition <- sapply(parts, function(x) x[2])
meta$time <- sapply(parts, function(x) x[3])

rownames(meta) <- meta$sample_name

# ---- 3. Filter for 0min samples ----
message("Filtering for 0min samples...")
keep_samples <- meta$time == "0min"
meta_sub <- meta[keep_samples, ]
counts_sub <- counts[, keep_samples]

meta_sub$condition <- factor(meta_sub$condition)
if ("WT" %in% levels(meta_sub$condition)) {
  meta_sub$condition <- relevel(meta_sub$condition, ref = "WT")
}

message(paste("Samples included:",
              paste(rownames(meta_sub), collapse = ", ")))

# ---- 4. Run DESeq2 ----
message("Running DESeq2...")
dds <- DESeqDataSetFromMatrix(countData = counts_sub,
                              colData = meta_sub,
                              design = ~ condition)

# Pre-filter: keep genes with >= 10 total counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds <- DESeq(dds)

res <- results(dds)
message("DESeq2 analysis complete.")

# ---- 5. Volcano plot ----
res_df <- as.data.frame(res)
res_df$gene_name <- rownames(res_df)

# Classify genes: padj < 0.05 and |log2FC| > 1
res_df$significant <- "No"
res_df$significant[which(res_df$padj < 0.05 &
                         res_df$log2FoldChange > 1)] <- "Up"
res_df$significant[which(res_df$padj < 0.05 &
                         res_df$log2FoldChange < -1)] <- "Down"

# Select top 10 up and top 10 down by lowest padj for labeling
top_up <- res_df[res_df$significant == "Up", ]
top_up <- top_up[order(top_up$padj), ][1:min(10, nrow(top_up)), ]

top_down <- res_df[res_df$significant == "Down", ]
top_down <- top_down[order(top_down$padj), ][1:min(10, nrow(top_down)), ]

label_genes <- rbind(top_up, top_down)
label_genes <- label_genes[!is.na(label_genes$gene_name), ]

p <- ggplot(res_df,
            aes(x = log2FoldChange, y = -log10(padj),
                color = significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Down" = "blue", "No" = "grey",
                                "Up" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",
             color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "black", alpha = 0.5) +
  theme_minimal() +
  labs(title = "Volcano Plot: WT vs EV (0min)",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text_repel(data = label_genes, aes(label = gene_name),
                  size = 3, box.padding = 0.5, point.padding = 0.5,
                  force = 2, max.overlaps = Inf,
                  show.legend = FALSE, color = "black")

# ---- 6. Save outputs ----
message("Saving results...")
ggsave("volcano_plot_0min.pdf", plot = p, width = 8, height = 6)
ggsave("volcano_plot_0min.png", plot = p, width = 8, height = 6)
write.csv(res_df, "deseq2_results_0min.csv", row.names = FALSE)

message("Done! Results saved in ", output_dir)
