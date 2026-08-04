# ============================================================
# Script:    run_go_enrichment_gradient.R
# Pipeline:  DESeq2 - GO Enrichment Analysis
# Purpose:   GO enrichment (BP, CC, MF) with gradient-colored bars
#            reflecting adjusted p-value (strong color = low p).
#            Terms sorted by -log10(p.adjust), top 10 per ontology.
# Input:     deseq2_results_0min.csv
# Output:    GO_Analysis_ClusterProfiler/GO_Enrichment_Gradient.pdf/.png
# Author:    Guchen-Weng
# ============================================================

suppressPackageStartupMessages({
  library(org.EcK12.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(cowplot)
  library(dplyr)
  library(scales)
})

# ---- Paths ----
input_file <- "deseq2_results_0min.csv"
output_dir <- "GO_Analysis_ClusterProfiler"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

setwd(output_dir)

# ---- 1. Read DESeq2 results ----
message("Reading DESeq2 results...")
res <- read.csv(input_file, stringsAsFactors = FALSE)

# ---- 2. Filter significant genes ----
sig_genes_df <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
message(paste("Number of significant genes:", nrow(sig_genes_df)))

if (nrow(sig_genes_df) == 0) {
  stop("No significant genes found with current thresholds.")
}

# ---- 3. Map gene symbols to Entrez IDs ----
message("Mapping Gene Symbols to Entrez IDs...")
gene_symbols <- sig_genes_df$gene_name

entrez_ids <- mapIds(org.EcK12.eg.db,
                     keys = gene_symbols,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")

valid_genes <- entrez_ids[!is.na(entrez_ids)]
message(paste("Mapped", length(valid_genes), "genes to Entrez IDs."))

# ---- 4. Run GO enrichment with gradient fill ----
ontologies <- c("BP", "CC", "MF")

# Low p-value (strong) color per ontology
ont_colors_low <- c("BP" = "#d74596", "CC" = "#95c94a", "MF" = "#984EA3")
# High p-value (light) color per ontology
ont_colors_high <- c("BP" = "#FADCEE", "CC" = "#EAF7D9", "MF" = "#EBD6F0")

plots_list <- list()

for (ont in ontologies) {
  message(paste("Running GO enrichment for:", ont))

  ego <- enrichGO(gene          = valid_genes,
                  OrgDb         = org.EcK12.eg.db,
                  keyType       = "ENTREZID",
                  ont           = ont,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2,
                  readable      = TRUE)

  if (is.null(ego) || nrow(ego) == 0) {
    message(paste("No enriched terms found for", ont))
    p <- ggplot() +
      theme_void() +
      labs(title = paste("No enriched terms for", ont))
  } else {
    res_df <- as.data.frame(ego)

    # Compute -log10(p.adjust) and select top 10
    res_df$log10P <- -log10(res_df$p.adjust)
    top_n <- 10
    res_df <- res_df[order(res_df$log10P, decreasing = TRUE), ]
    if (nrow(res_df) > top_n) {
      res_df <- res_df[1:top_n, ]
    }

    # Reorder so highest log10P is at the top of the plot
    res_df <- res_df[order(res_df$log10P, decreasing = FALSE), ]
    res_df$Description <- factor(res_df$Description,
                                 levels = res_df$Description)

    p <- ggplot(res_df,
                aes(x = log10P, y = Description, fill = p.adjust)) +
      geom_col(width = 0.7) +
      scale_fill_gradient(low = ont_colors_low[ont],
                          high = ont_colors_high[ont],
                          name = "Adj. P-value",
                          guide = guide_colorbar(reverse = TRUE)) +
      theme_minimal() +
      labs(title = paste(ont),
           x = "-log10(Adjusted P-value)",
           y = NULL) +
      theme(axis.text.y = element_text(size = 10, color = "black"),
            axis.text.x = element_text(size = 10, color = "black"),
            plot.title = element_text(size = 12, face = "bold",
                                      hjust = -0.1),
            legend.position = "right")
  }

  plots_list[[ont]] <- p
}

# ---- 5. Combine plots ----
message("Combining plots...")
combined_plot <- plot_grid(plots_list[["BP"]],
                           plots_list[["CC"]],
                           plots_list[["MF"]],
                           ncol = 1,
                           align = "v",
                           rel_heights = c(1, 1, 1))

# ---- 6. Save ----
output_pdf <- "GO_Enrichment_Gradient.pdf"
output_png <- "GO_Enrichment_Gradient.png"

ggsave(output_pdf, combined_plot, width = 10, height = 12)
ggsave(output_png, combined_plot, width = 10, height = 12)

message(paste("Saved combined plot to", output_pdf))
