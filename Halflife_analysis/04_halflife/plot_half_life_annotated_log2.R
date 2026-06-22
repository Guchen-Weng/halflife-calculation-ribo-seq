# ============================================================
# Script:    plot_half_life_annotated_log2.R
# Pipeline:  Step 04 - Half-life Visualization
# Purpose:   Violin+boxplot of log2 half-life distributions by
#            Condition (EV/WT) and motif presence (True/False),
#            annotated with Wilcoxon rank-sum test p-values
# Input:     ../../half_decay/half_life_results.csv
# Output:    ../output/half_life_violin_annotated_log2.pdf
# Author:    Guchen-Weng
# Date:      2025-01
# ============================================================

library(ggplot2)
library(grid)

# ---- Paths ----
input_file <- "../../half_decay/half_life_results.csv"
output_dir <- "../output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output_plot <- file.path(output_dir, "half_life_violin_annotated_log2.pdf")

# ---- Read and filter data ----
data <- read.csv(input_file, stringsAsFactors = FALSE)
data <- data[!is.na(data$half_life), ]
data <- data[data$half_life <= 2880, ]
data <- data[data$half_life > 0, ]

# Log2 transformation
data$log2_half_life <- log2(data$half_life)
data <- data[data$log2_half_life >= -1, ]

# Factor encoding
data$has_motif <- factor(data$has_motif,
                         levels = c("FALSE", "TRUE"),
                         labels = c("False", "True"))
data$Condition <- factor(data$Condition, levels = c("EV", "WT"))

# Create combined Group factor
data$Group <- paste0(data$Condition, "_", data$has_motif)
data$Group <- factor(data$Group,
                     levels = c("EV_False", "EV_True", "WT_False", "WT_True"))

# ---- Counts for legend ----
counts <- aggregate(log2_half_life ~ Condition + has_motif,
                    data = data, length)
colnames(counts)[3] <- "count"
data <- merge(data, counts, by = c("Condition", "has_motif"))

# Legend labels
legend_map <- unique(data[, c("Group", "Condition", "has_motif", "count")])
legend_map$Label <- paste0(legend_map$Condition, ", Motif: ",
                           legend_map$has_motif,
                           " (n=", legend_map$count, ")")
group_labels <- setNames(legend_map$Label, legend_map$Group)

# ---- Custom colors ----
my_colors <- c("WT_True"  = "#d74596",
               "WT_False" = "#EBA2CA",
               "EV_True"  = "#95c94a",
               "EV_False" = "#CAE4A4")

# ---- Compute Wilcoxon p-values for annotation ----
annotations <- data.frame()
max_y <- max(data$log2_half_life, na.rm = TRUE)
y_pos <- max_y * 1.05

for (cond in c("EV", "WT")) {
  sub_data <- data[data$Condition == cond, ]
  res <- wilcox.test(log2_half_life ~ has_motif, data = sub_data)
  p_val <- res$p.value
  p_label <- paste0("p = ", formatC(p_val, format = "e", digits = 2))

  annotations <- rbind(annotations, data.frame(
    Condition = cond,
    log2_half_life = y_pos,
    label = p_label
  ))
}
annotations$Condition <- factor(annotations$Condition, levels = c("EV", "WT"))

# ---- Plot ----
p <- ggplot(data, aes(x = Condition, y = log2_half_life, fill = Group)) +
  geom_violin(width = 0.7, position = position_dodge(width = 0.9),
              trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.9),
               outlier.shape = NA, alpha = 0.5, show.legend = FALSE) +

  # p-value annotations
  geom_text(data = annotations,
            aes(x = Condition, y = log2_half_life, label = label),
            inherit.aes = FALSE, size = 5, fontface = "bold") +

  labs(y = "Log2 Half-life", x = "Condition",
       fill = "Group (Count)", title = "Log2 Half-life Distribution") +
  theme_bw() +
  theme(
    legend.position = "right",
    text = element_text(size = 12)
  ) +

  scale_fill_manual(values = my_colors, labels = group_labels) +
  guides(fill = guide_legend(override.aes = list(color = "black"))) +
  scale_y_continuous(breaks = seq(-4, 20, by = 2)) +
  coord_cartesian(ylim = c(-3, 20))

# ---- Save ----
ggsave(output_plot, plot = p, width = 10, height = 8)
print(paste("Plot saved to", output_plot))
