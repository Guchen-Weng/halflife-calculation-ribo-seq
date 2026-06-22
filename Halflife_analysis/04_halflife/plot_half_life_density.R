# ============================================================
# Script:    plot_half_life_density.R
# Pipeline:  Step 04 - Half-life Visualization
# Purpose:   Density plot of half-life values by Condition + motif,
#            with peak annotations (vertical dashed lines + labels)
# Input:     ../../half_decay/half_life_results.csv
# Output:    ../output/half_life_density.pdf
# Author:    Guchen-Weng
# Date:      2025-01
# ============================================================

library(ggplot2)

# ---- Paths ----
input_file <- "../../half_decay/half_life_results.csv"
output_dir <- "../output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output_plot <- file.path(output_dir, "half_life_density.pdf")

# ---- Read and filter data ----
data <- read.csv(input_file, stringsAsFactors = FALSE)
data <- data[!is.na(data$half_life), ]
data <- data[data$half_life <= 2880, ]

# Factor encoding
data$has_motif <- factor(data$has_motif,
                         levels = c("FALSE", "TRUE"),
                         labels = c("False", "True"))
data$Condition <- factor(data$Condition, levels = c("EV", "WT"))

# Create combined Group factor
data$Group <- paste0(data$Condition, "_", data$has_motif)
data$Group <- factor(data$Group,
                     levels = c("EV_False", "EV_True", "WT_False", "WT_True"))

# ---- Find density peaks per group ----
get_density_peak <- function(x) {
  if (length(x) < 2) return(data.frame(x = NA, y = NA))
  d <- density(x, from = min(x), to = max(x))
  max_y <- max(d$y)
  max_x <- d$x[which.max(d$y)]
  return(data.frame(x = max_x, y = max_y))
}

peaks_list <- lapply(split(data, data$Group), function(df) {
  peak <- get_density_peak(df$half_life)
  peak$Group <- unique(df$Group)
  return(peak)
})
peaks <- do.call(rbind, peaks_list)
rownames(peaks) <- NULL

# ---- Custom colors ----
my_colors <- c("WT_True"  = "#d74596",
               "WT_False" = "#EBA2CA",
               "EV_True"  = "#00552E",
               "EV_False" = "#95c94a")

# ---- Plot ----
p <- ggplot(data, aes(x = half_life, color = Group, fill = Group)) +
  geom_density(alpha = 0.3) +

  # Vertical dashed lines from peak to x-axis
  geom_segment(data = peaks,
               aes(x = x, xend = x, y = 0, yend = y, color = Group),
               linetype = "dashed", size = 0.5) +

  # Peak value labels
  geom_text(data = peaks,
            aes(x = x, y = y, label = round(x, 1), color = Group),
            vjust = -0.5, fontface = "bold", show.legend = FALSE) +

  labs(x = "Half-life (min)", y = "Density",
       title = "Half-life Density Distribution") +
  theme_bw() +
  theme(
    legend.position = "right",
    text = element_text(size = 12)
  ) +

  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +

  coord_cartesian(xlim = c(0, 500))

# ---- Save ----
ggsave(output_plot, plot = p, width = 10, height = 6)
print(paste("Plot saved to", output_plot))
