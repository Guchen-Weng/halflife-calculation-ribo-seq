# ============================================================
# Script:    run_go_enrichment_combined.R
# Pipeline:  DESeq2 - GO Enrichment Analysis
# Purpose:   GO enrichment (BP, CC, MF) using clusterProfiler on
#            significant genes from DESeq2 results. Produces a
#            combined barplot of top 10 terms per ontology, stacked
#            vertically with fixed per-ontology colors.
# Input:     deseq2_results_0min.csv
# Output:    GO_Analysis_ClusterProfiler/GO_Enrichment_Combined.pdf/.png
# Author:    Guchen-Weng
# ============================================================

suppressPackageStartupMessages({
  library(org.EcK12.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(cowplot)
  library(dplyr)
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

# ---- 2. Filter significant genes (padj < 0.05 & |log2FC| > 1) ----
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

# ---- 4. Run GO enrichment (BP, CC, MF) ----
ontologies <- c("BP", "CC", "MF")
ont_colors <- c("BP" = "#d74596", "CC" = "#95c94a", "MF" = "#984EA3")

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
    # Extract top 10 terms
    res_df <- as.data.frame(ego)
    top_n <- 10
    if (nrow(res_df) > top_n) {
      res_df <- res_df[1:top_n, ]
    }

    # Order by Count (descending) -> factor levels reversed for ggplot
    res_df$Description <- factor(res_df$Description,
                                 levels = rev(res_df$Description))

    p <- ggplot(res_df, aes(x = Count, y = Description)) +
      geom_bar(stat = "identity", fill = ont_colors[ont], width = 0.7) +
      theme_minimal() +
      labs(title = paste("GO Enrichment:", ont),
           x = "Gene Count",
           y = NULL) +
      theme(axis.text.y = element_text(size = 10),
            plot.title = element_text(size = 12, face = "bold"))
  }

  plots_list[[ont]] <- p
}

# ---- 5. Combine plots (stacked vertically) ----
message("Combining plots...")
combined_plot <- plot_grid(plots_list[["BP"]],
                           plots_list[["CC"]],
                           plots_list[["MF"]],
                           ncol = 1,
                           align = "v",
                           labels = "AUTO")

# ---- 6. Save ----
output_pdf <- "GO_Enrichment_Combined.pdf"
output_png <- "GO_Enrichment_Combined.png"

ggsave(output_pdf, combined_plot, width = 10, height = 15)
ggsave(output_png, combined_plot, width = 10, height = 15)

message(paste("Saved combined plot to", output_pdf))
