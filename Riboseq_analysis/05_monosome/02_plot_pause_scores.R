#!/usr/bin/env Rscript
# ============================================================
# Script:    02_plot_pause_scores.R
# Pipeline:  Monosome Pause Score Visualization
# Purpose:   Plot EV vs WT codon pause scores (A/P/E sites) from
#            the output of 01_pause_score_calculation.R.
#            Dots = group mean, error bars = SE across replicates.
# Input:     PauseScore_fixedStrand_len26-32/
#              *_pause_scores_codon_A_P_E.csv
# Output:    codon_pause_EV_vs_WT_APE_minRPN0.1.png/.pdf
#            codon_pause_EV_vs_WT_APE_minRPN0.1_summary.csv
# Author:    Guchen-Weng
# Date:      2025-01
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(stringr)
})

# =====================================================================
# 0) Input files
# =====================================================================
workdir <- "PauseScore_fixedStrand_len26-32"
if (!dir.exists(workdir)) stop("Input directory not found: ", workdir)
setwd(workdir)

files <- c(
  "EV_1_cds_annotation_with_PEA_offset15_fixedStrand_pause_scores_codon_A_P_E.csv",
  "EV_2_cds_annotation_with_PEA_offset15_fixedStrand_pause_scores_codon_A_P_E.csv",
  "WT_1_cds_annotation_with_PEA_offset15_fixedStrand_pause_scores_codon_A_P_E.csv",
  "WT_2_cds_annotation_with_PEA_offset15_fixedStrand_pause_scores_codon_A_P_E.csv"
)
stopifnot(all(file.exists(files)))

# =====================================================================
# 1) Parameters
# =====================================================================
minRPN <- 0.1
out_prefix <- paste0("codon_pause_EV_vs_WT_APE_minRPN", minRPN)

# Fixed codon order: AAA, AAC, AAG, AAT, ACA, ... , TTT
codon_levels <- as.vector(outer(c("A","C","G","T"),
  outer(c("A","C","G","T"), c("A","C","G","T"), paste0), paste0))

# =====================================================================
# 2) Read and annotate each table
# =====================================================================
read_one <- function(f) {
  dt <- fread(f)

  if (!("site" %in% names(dt))) stop("Missing column 'site' in: ", f)
  if (!("codon" %in% names(dt))) {
    if ("codon_nt" %in% names(dt)) setnames(dt, "codon_nt", "codon")
    else stop("Missing column 'codon' or 'codon_nt' in: ", f)
  }
  if (!("pause_mean" %in% names(dt))) {
    if ("Mean_pause" %in% names(dt)) setnames(dt, "Mean_pause", "pause_mean")
    else stop("Missing column 'pause_mean' or 'Mean_pause' in: ", f)
  }
  if (!("n" %in% names(dt))) {
    if ("N_instances" %in% names(dt)) setnames(dt, "N_instances", "n")
    else dt[, n := NA_integer_]
  }

  samp <- str_replace(basename(f), "_pause_scores_codon_A_P_E\\.csv$", "")
  grp  <- ifelse(str_detect(samp, "^EV"), "EV", "WT")
  dt[, `:=`(sample = samp, group = grp)]
  dt
}

dt_all <- rbindlist(lapply(files, read_one), use.names = TRUE, fill = TRUE)

# =====================================================================
# 3) Optional minRPN filtering
# =====================================================================
rpn_col <- intersect(names(dt_all),
  c("reads_per_nt", "readpernt", "readspernt", "read_per_nt"))
if (length(rpn_col) > 0) {
  rpn_col <- rpn_col[1]
  dt_all <- dt_all[get(rpn_col) >= minRPN]
  message("Applied minRPN filter using column: ", rpn_col)
} else {
  message("No RPN column found -> skip minRPN filtering.")
}

# =====================================================================
# 4) Standardize factor levels (A/P/E sites only)
# =====================================================================
dt_all <- dt_all[site %in% c("A", "P", "E")]
dt_all[, codon := toupper(codon)]
dt_all <- dt_all[nchar(codon) == 3 & !grepl("N", codon)]
dt_all[, codon := factor(codon, levels = codon_levels)]
dt_all[, site  := factor(site, levels = c("A", "E", "P"))]
dt_all[, group := factor(group, levels = c("EV", "WT"))]

# =====================================================================
# 5) Mean +/- SE across replicates per codon/site/group
# =====================================================================
sum_dt <- dt_all[, .(
  mean_pause = mean(pause_mean, na.rm = TRUE),
  sd_pause   = sd(pause_mean,  na.rm = TRUE),
  n_rep      = sum(!is.na(pause_mean)),
  se_pause   = sd(pause_mean, na.rm = TRUE) /
               sqrt(sum(!is.na(pause_mean)))
), by = .(site, codon, group)]

sum_dt <- sum_dt[!is.na(codon)]

# =====================================================================
# 6) Plot: group means with SE error bars, faceted by A/E/P
# =====================================================================
p <- ggplot(sum_dt, aes(x = codon, y = mean_pause, color = group)) +
  geom_errorbar(aes(ymin = mean_pause - se_pause,
                    ymax = mean_pause + se_pause),
                width = 0.25, linewidth = 0.4,
                position = position_dodge(width = 0.5)) +
  geom_point(size = 2.0, position = position_dodge(width = 0.5)) +
  facet_wrap(~ site, ncol = 1, scales = "free_y") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_line(linewidth = 0.2),
    panel.grid.minor.x = element_blank()
  ) +
  labs(
    title = paste0("Codon pause EV vs WT (A/P/E sites) | minRPN ", minRPN),
    subtitle = "Dots = group mean; error bars = SE (EV n=2, WT n=2)",
    x = "Codon (AAA to TTT)",
    y = "Pause score (mean +/- SE)",
    color = "Group"
  )

print(p)

# =====================================================================
# 7) Save outputs
# =====================================================================
ggsave(paste0(out_prefix, ".png"), p, width = 16, height = 6.5, dpi = 300)
ggsave(paste0(out_prefix, ".pdf"), p, width = 16, height = 6.5)
fwrite(sum_dt, paste0(out_prefix, "_summary.csv"))

message("Done. Wrote: ", out_prefix, ".png/.pdf and summary csv")
