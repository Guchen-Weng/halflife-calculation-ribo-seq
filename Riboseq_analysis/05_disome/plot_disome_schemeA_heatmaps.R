#!/usr/bin/env Rscript
# ============================================================
# Script:    plot_disome_schemeA_heatmaps.R
# Pipeline:  Disome Analysis - Visualization
# Purpose:   Generate "Scheme A" plots from the disome cache RDS:
#            top panel = mean dZ curve with SEM ribbon,
#            bottom panel = dZ heatmap by read length.
#            Produces per-target PNG/PDF and a combined multi-page
#            PDF for all targets (stop + motifs).
# Input:     cache_stop_and_motifs_hm_dz.rds
# Output:    SchemeA_avg_plus_heatmap_36_60nt/
#              *_schemeA_avg_plus_heatmap_36_60nt.pdf / .png
#              ALL_targets_average_dZ_36_60nt.csv
# Author:    Guchen-Weng
# Date:      2025-01
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# =====================================================================
# 0) Path configuration
# =====================================================================
cache_rds <- "cache_stop_and_motifs_hm_dz.rds"

len_min <- 36
len_max <- 60
out_dir <- paste0("SchemeA_avg_plus_heatmap_", len_min, "_", len_max, "nt")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =====================================================================
# 1) Read the RDS cache
# =====================================================================
out <- readRDS(cache_rds)

# =====================================================================
# 2) Auto-detect available targets
# =====================================================================
targets <- c()

if ("stop" %in% names(out))  targets <- c(targets, "stop")
if ("stops" %in% names(out)) targets <- c(targets, "stops")

if ("motifs" %in% names(out)) {
  targets <- c(targets, names(out$motifs))
}

targets <- unique(targets)

cat("Targets found in cache:\n")
print(targets)

preferred_order <- c("stop", "stops", "ATA", "TAA", "TAAG", "CTC", "GCT")
targets <- preferred_order[preferred_order %in% targets]

cat("Targets to plot:\n")
print(targets)

# =====================================================================
# 3) Helper: retrieve object by target name
# =====================================================================
get_target_obj <- function(out, target_name) {
  if (target_name == "stop" && "stop" %in% names(out)) {
    return(out$stop)
  }
  if (target_name == "stops" && "stops" %in% names(out)) {
    return(out$stops)
  }
  if ("motifs" %in% names(out) && target_name %in% names(out$motifs)) {
    return(out$motifs[[target_name]])
  }
  stop("Cannot find target: ", target_name)
}

# =====================================================================
# 4) Plot-saving helpers
# =====================================================================
save_plot_pdf <- function(plot_obj, filename, width = 10.5, height = 7.2) {
  grDevices::pdf(filename, width = width, height = height, useDingbats = FALSE)
  print(plot_obj)
  grDevices::dev.off()
}

save_plot_png <- function(plot_obj, filename, width = 10.5, height = 7.2, res = 300) {
  grDevices::png(filename, width = width, height = height,
                 units = "in", res = res)
  print(plot_obj)
  grDevices::dev.off()
}

# =====================================================================
# 5) Generate a Scheme A plot for a single target
#    Top panel: mean dZ curve with SEM ribbon
#    Bottom panel: heatmap
# =====================================================================
make_schemeA_plot <- function(out, target_name, len_min = 36, len_max = 60) {

  obj <- get_target_obj(out, target_name)

  if (!("dz5" %in% names(obj)) || !("dz3" %in% names(obj))) {
    stop("Target ", target_name, " does not contain dz5 and dz3")
  }

  dz5 <- as.data.table(obj$dz5)
  dz3 <- as.data.table(obj$dz3)

  needed_cols <- c("read_length", "dist", "dz")
  if (!all(needed_cols %in% names(dz5))) {
    stop("dz5 of ", target_name, " missing columns: ",
         paste(setdiff(needed_cols, names(dz5)), collapse = ", "))
  }
  if (!all(needed_cols %in% names(dz3))) {
    stop("dz3 of ", target_name, " missing columns: ",
         paste(setdiff(needed_cols, names(dz3)), collapse = ", "))
  }

  dz5 <- dz5[read_length >= len_min & read_length <= len_max]
  dz3 <- dz3[read_length >= len_min & read_length <= len_max]

  dz5[, end_type := "5' end"]
  dz3[, end_type := "3' end"]

  # ---- Average curve (mean dZ +/- SEM) ----
  avg5 <- dz5[, .(
    mean_dz = mean(dz, na.rm = TRUE),
    sem_dz  = sd(dz, na.rm = TRUE) / sqrt(sum(!is.na(dz)))
  ), by = .(dist)]
  avg5[, end_type := "5' end"]

  avg3 <- dz3[, .(
    mean_dz = mean(dz, na.rm = TRUE),
    sem_dz  = sd(dz, na.rm = TRUE) / sqrt(sum(!is.na(dz)))
  ), by = .(dist)]
  avg3[, end_type := "3' end"]

  avg_dt <- rbindlist(list(avg3, avg5), use.names = TRUE)

  # ---- Heatmap data ----
  hm_dt <- rbindlist(list(
    dz3[, .(read_length, dist, dz, end_type)],
    dz5[, .(read_length, dist, dz, end_type)]
  ), use.names = TRUE)

  hm_dt[, read_length := as.integer(read_length)]

  # ---- Top panel: mean dZ curve ----
  p_avg <- ggplot(avg_dt, aes(x = dist, y = mean_dz)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey70", linewidth = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey35", linewidth = 0.6) +
    geom_ribbon(aes(ymin = mean_dz - sem_dz, ymax = mean_dz + sem_dz),
                fill = "#E91E63", alpha = 0.18, color = NA) +
    geom_line(color = "#E91E63", linewidth = 1) +
    facet_wrap(~ end_type, scales = "free_x", nrow = 1) +
    labs(
      title = paste0("Average dZ signal around ", target_name, " sites"),
      subtitle = paste0("dZ = WT - EV | Read lengths ", len_min, "-", len_max, " nt"),
      x = NULL, y = "Mean dZ"
    ) +
    theme_classic(base_size = 15) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      strip.background = element_rect(fill = "white", color = "black",
                                      linewidth = 1),
      strip.text = element_text(face = "bold", size = 13),
      axis.title.y = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8)
    )

  # ---- Bottom panel: heatmap ----
  p_hm <- ggplot(hm_dt, aes(x = dist, y = read_length, fill = dz)) +
    geom_tile() +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "black", linewidth = 0.5) +
    facet_wrap(~ end_type, scales = "free_x", nrow = 1) +
    scale_fill_gradient2(
      low = "#3B4CC0", mid = "white", high = "#B40426",
      midpoint = 0, limits = c(-2, 2), oob = scales::squish,
      name = "dZ\n(WT-EV)"
    ) +
    scale_y_continuous(
      breaks = seq(len_min, len_max, by = 4), expand = c(0, 0)
    ) +
    labs(
      x = paste0("Distance from ", target_name, " anchor (nt)"),
      y = "Read length (nt)"
    ) +
    theme_classic(base_size = 15) +
    theme(
      strip.background = element_rect(fill = "white", color = "black",
                                      linewidth = 1),
      strip.text = element_text(face = "bold", size = 13),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      legend.position = "right"
    )

  # ---- Combine panels ----
  p_all <- p_avg / p_hm + plot_layout(heights = c(1, 1.2))

  return(list(plot = p_all, avg_dt = avg_dt, hm_dt = hm_dt))
}

# =====================================================================
# 6) Generate and save plots for all targets
# =====================================================================
plot_list <- list()
avg_all_list <- list()

for (tg in targets) {
  cat("Plotting:", tg, "\n")

  res <- make_schemeA_plot(out, tg, len_min = len_min, len_max = len_max)

  plot_list[[tg]] <- res$plot

  avg_tmp <- copy(res$avg_dt)
  avg_tmp[, target := tg]
  avg_all_list[[tg]] <- avg_tmp

  prefix <- file.path(out_dir,
    paste0(tg, "_schemeA_avg_plus_heatmap_", len_min, "_", len_max, "nt"))

  save_plot_pdf(res$plot, paste0(prefix, ".pdf"),
                width = 10.5, height = 7.2)
  save_plot_png(res$plot, paste0(prefix, ".png"),
                width = 10.5, height = 7.2, res = 300)
}

# =====================================================================
# 7) Multi-page PDF (all targets)
# =====================================================================
multi_pdf <- file.path(out_dir,
  paste0("ALL_targets_schemeA_avg_plus_heatmap_", len_min, "_", len_max, "nt.pdf"))
grDevices::pdf(multi_pdf, width = 10.5, height = 7.2, useDingbats = FALSE)

for (tg in names(plot_list)) {
  print(plot_list[[tg]])
}

grDevices::dev.off()

# =====================================================================
# 8) Export mean-curve data
# =====================================================================
avg_all <- rbindlist(avg_all_list, use.names = TRUE, fill = TRUE)

fwrite(avg_all,
  file = file.path(out_dir,
    paste0("ALL_targets_average_dZ_", len_min, "_", len_max, "nt.csv")))

cat("\nDone!\n")
cat("Output directory:\n", out_dir, "\n")
cat("Multi-page PDF:\n", multi_pdf, "\n")
