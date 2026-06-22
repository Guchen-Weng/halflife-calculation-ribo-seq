# ============================================================
# Script:    plot_half_life_ecdf.R
# Pipeline:  Step 04 - Half-life Visualization
# Purpose:   Empirical cumulative distribution function (ECDF) plot
#            comparing EV vs. WT half-life distributions
# Input:     ../../half_decay/half_life_results.csv
# Output:    ../output/half_life_ecdf.pdf
# Author:    Guchen-Weng
# Date:      2025-01
# ============================================================

library(ggplot2)

# ---- Paths ----
input_file <- "../../half_decay/half_life_results.csv"
output_dir <- "../output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output_plot <- file.path(output_dir, "half_life_ecdf.pdf")

# ---- Read and filter data ----
data <- read.csv(input_file, stringsAsFactors = FALSE)
data <- data[!is.na(data$half_life), ]
data <- data[data$half_life <= 2880, ]

# Factor encoding
data$Condition <- factor(data$Condition, levels = c("EV", "WT"))

# ---- Custom colors ----
my_colors <- c("WT" = "#d74596", "EV" = "#00552E")

# ---- Plot ----
p <- ggplot(data, aes(x = half_life, color = Condition)) +
  stat_ecdf(geom = "step", linewidth = 1) +

  labs(x = "Half-life (min)", y = "Cumulative Probability",
       title = "Cumulative Distribution of Half-lives") +
  theme_bw() +
  theme(
    legend.position = "right",
    text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5)
  ) +

  scale_color_manual(values = my_colors) +
  coord_cartesian(xlim = c(0, 100))

# ---- Save ----
ggsave(output_plot, plot = p, width = 8, height = 6)
print(paste("Plot saved to", output_plot))
