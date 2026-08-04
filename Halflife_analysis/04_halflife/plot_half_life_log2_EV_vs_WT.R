# ============================================================
# Script:    plot_half_life_log2_EV_vs_WT.R
# Pipeline:  Step 04 - Half-life Visualization
# Purpose:   Violin+boxplot comparing log2 half-life between EV and
#            WT (ignoring motif status), with sample counts in legend
# Input:     ../../half_decay/half_life_results.csv
# Output:    ../output/half_life_violin_log2_EV_vs_WT.pdf
# Author:    Guchen-Weng
# ============================================================

library(ggplot2)
library(grid)

# ---- Paths ----
input_file <- "../../half_decay/half_life_results.csv"
output_dir <- "../output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output_plot <- file.path(output_dir, "half_life_violin_log2_EV_vs_WT.pdf")

# ---- Read and filter data ----
data <- read.csv(input_file, stringsAsFactors = FALSE)
data <- data[!is.na(data$half_life), ]
data <- data[data$half_life <= 2880, ]
data <- data[data$half_life > 0, ]

# Log2 transformation
data$log2_half_life <- log2(data$half_life)
data <- data[data$log2_half_life >= -1, ]

# Factor encoding
data$Condition <- factor(data$Condition, levels = c("EV", "WT"))

# ---- Counts for legend ----
counts <- aggregate(log2_half_life ~ Condition, data = data, length)
colnames(counts)[2] <- "count"
data <- merge(data, counts, by = "Condition")

legend_map <- unique(data[, c("Condition", "count")])
legend_map$Label <- paste0(legend_map$Condition,
                           " (n=", legend_map$count, ")")
group_labels <- setNames(legend_map$Label, legend_map$Condition)

# ---- Custom colors ----
my_colors <- c("WT" = "#d74596", "EV" = "#95c94a")

# ---- Plot ----
p <- ggplot(data, aes(x = Condition, y = log2_half_life, fill = Condition)) +
  geom_violin(width = 0.7, trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, outlier.shape = NA,
               alpha = 0.5, show.legend = FALSE) +

  labs(y = "Log2 Half-life", x = "Condition",
       fill = "Condition (Count)",
       title = "Log2 Half-life Distribution (EV vs WT)") +
  theme_bw() +
  theme(
    legend.position = "right",
    text = element_text(size = 12),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black")
  ) +

  scale_y_continuous(breaks = seq(0, 20, by = 2)) +
  coord_cartesian(ylim = c(-3, 20)) +
  scale_fill_manual(values = my_colors, labels = group_labels) +
  guides(fill = guide_legend(override.aes = list(color = "black")))

# ---- Save ----
ggsave(output_plot, plot = p, width = 8, height = 8)
print(paste("Plot saved to", output_plot))
