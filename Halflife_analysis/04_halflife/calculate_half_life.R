# ============================================================
# Script:    calculate_half_life.R
# Pipeline:  Step 04 - Half-life Calculation
# Purpose:   Estimate mRNA half-life per gene per condition (WT/EV)
#            by fitting exponential decay curves using NLS, with
#            LM fallback on log-transformed data when NLS fails.
#            half_life = ln(2) / k
# Input:     ../../UCSC_db/gene_tpm_with_motifs.csv
# Output:    ../../half_decay/half_life_results.csv
# Author:    Guchen-Weng
# ============================================================

# ---- Paths ----
# Input: gene TPM matrix with motif annotations
# Expected columns: Geneid, chromosome, start, end, strand, biotype,
#   description, gene_name, motif_count_UAAAG, motif_count_UAAGA,
#   total_motifs, has_motif, plus sample columns
input_file <- "../../UCSC_db/gene_tpm_with_motifs.csv"

# Output: half-life estimates per gene per condition
output_dir <- "../../half_decay"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output_file <- file.path(output_dir, "half_life_results.csv")

# ---- Read data ----
print(paste("Reading data from:", input_file))
data <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE)

# ---- fit_decay function ----
# Fits exponential decay: expression ~ a * exp(-k * time)
# method = "nls": non-linear least squares (preferred)
# method = "lm":  linear model on log-transformed data (fallback)
fit_decay <- function(dt = NULL,
                      time_col = "time",
                      expression_col = "expression",
                      method = "nls",
                      verbose = FALSE) {
  if (!time_col %in% names(dt))
    stop(paste("Time column", time_col, "not found in table!"))
  if (!expression_col %in% names(dt))
    stop(paste("Expression column", expression_col, "not found in table!"))
  if (!is.numeric(dt[[time_col]]))
    stop(paste(time_col, "must be numeric!"))
  if (!is.numeric(dt[[expression_col]]))
    stop(paste(expression_col, "must be numeric!"))

  if (method == "nls") {
    fml = as.formula(paste(expression_col, "~ I(a * exp(-k *", time_col, "))"))

    # Filter invalid values
    valid_data <- dt[is.finite(dt[[expression_col]]) &
                     !is.na(dt[[expression_col]]), ]
    if (nrow(valid_data) == 0) max_expr <- 0
    else max_expr <- max(valid_data[[expression_col]])

    if (nrow(valid_data) < 3) {
        fit = NA
        a = as.numeric(NA)
        a_se = as.numeric(NA)
        k = as.numeric(NA)
        k_se = as.numeric(NA)
    } else if (inherits(try(nls(
      formula =  fml,
      start = list(a = max_expr, k = 0.1),
      data = dt,
      control = nls.control(maxiter = 100, warnOnly = TRUE)
    ),
    silent = !verbose)
    , "try-error")) {
      fit = NA
      a = as.numeric(NA)
      a_se = as.numeric(NA)
      k = as.numeric(NA)
      k_se = as.numeric(NA)
    } else {
      fit = nls(formula =  fml,
                start = list(a = max_expr, k = 0.1),
                data = dt,
                control = nls.control(maxiter = 100, warnOnly = TRUE))

      smr_try <- try(summary(fit), silent = !verbose)

      if (inherits(smr_try, "try-error")) {
         fit = NA
         a = as.numeric(NA)
         a_se = as.numeric(NA)
         k = as.numeric(NA)
         k_se = as.numeric(NA)
      } else {
         smr = smr_try
         a = smr$parameters["a", "Estimate"]
         a_se = smr$parameters["a", "Std. Error"]
         k = smr$parameters["k", "Estimate"]
         k_se = smr$parameters["k", "Std. Error"]
      }
    }
  } else if (method == "lm") {
    fml = as.formula(paste("log(", expression_col, ")",  "~", time_col))

    # Filter out <= 0 values (log undefined)
    dt_lm <- dt[!is.na(dt[[expression_col]]) &
                is.finite(dt[[expression_col]]) &
                dt[[expression_col]] > 0, ]

    if (nrow(dt_lm) < 2) {
        fit = NA
        a = as.numeric(NA)
        a_se = as.numeric(NA)
        k = as.numeric(NA)
        k_se = as.numeric(NA)
    } else if (inherits(try(lm(
      formula =  fml,
      data = dt_lm,
      na.action = na.exclude
    ),
    silent = !verbose)
    , "try-error")) {
      fit = NA
      a = as.numeric(NA)
      a_se = as.numeric(NA)
      k = as.numeric(NA)
      k_se = as.numeric(NA)
    } else {
      fit = lm(formula =  fml,
               data = dt_lm,
               na.action = na.exclude)
      smr = summary(fit)

      intercept_est <- try(smr$coefficients["(Intercept)", "Estimate"],
                           silent = !verbose)
      intercept_se  <- try(smr$coefficients["(Intercept)", "Std. Error"],
                           silent = !verbose)
      slope_est     <- try(smr$coefficients[time_col, "Estimate"],
                           silent = !verbose)
      slope_se      <- try(smr$coefficients[time_col, "Std. Error"],
                           silent = !verbose)

      a    = if (inherits(intercept_est, "try-error")) as.numeric(NA)
             else exp(intercept_est)  # intercept = log(a)
      a_se = if (inherits(intercept_se, "try-error")) as.numeric(NA)
             else intercept_se
      k    = if (inherits(slope_est, "try-error")) as.numeric(NA)
             else -slope_est          # slope = -k
      k_se = if (inherits(slope_se, "try-error")) as.numeric(NA)
             else slope_se
    }
  } else {
    stop("Undefined method!")
  }

  coef <- data.frame(
      a = a,
      a_se = a_se,
      k = k,
      k_se = k_se,
      method = method,
      stringsAsFactors = FALSE)

  return(list(
    fit = fit,
    coef = coef))
}

# ---- Identify sample columns ----
non_sample_cols <- c("Geneid", "chromosome", "start", "end", "strand",
                     "biotype", "description", "gene_name",
                     "motif_count_UAAAG", "motif_count_UAAGA",
                     "total_motifs", "has_motif")
sample_cols <- setdiff(colnames(data), non_sample_cols)

# ---- Parse sample metadata from column names ----
# Expected format: Replicate_Condition_Time_ID (e.g., A1_WT_0min_1)
sample_info <- data.frame(Sample = sample_cols, stringsAsFactors = FALSE)

# Extract Condition (WT/EV)
sample_info$Condition <- ifelse(grepl("_WT_", sample_cols), "WT",
                                ifelse(grepl("_EV_", sample_cols), "EV", NA))

# Extract Time (0min, 2min, 6min, 12min)
sample_info$TimeStr <- regmatches(sample_cols, regexpr("[0-9]+min", sample_cols))
sample_info$Time <- as.numeric(gsub("min", "", sample_info$TimeStr))

# Remove unrecognized columns
sample_info <- sample_info[!is.na(sample_info$Condition) &
                           !is.na(sample_info$Time), ]

print(paste("Identified", nrow(sample_info), "valid samples."))
print(table(sample_info$Condition, sample_info$Time))

# ---- Fit half-life for each gene ----
print("Calculating half-life for genes...")
genes <- data$Geneid
n_genes <- length(genes)

process_gene <- function(i) {
    gene_id <- data$Geneid[i]
    row_vals <- data[i, sample_info$Sample]

    # Build per-gene data frame
    df_gene <- data.frame(
        Sample = names(row_vals),
        expression = as.numeric(row_vals),
        stringsAsFactors = FALSE
    )

    df_gene <- merge(df_gene, sample_info, by = "Sample")

    res_rows <- list()

    # Fit separately for WT and EV
    for (cond in c("WT", "EV")) {
        sub_df <- df_gene[df_gene$Condition == cond, ]

        # Try NLS first
        res <- fit_decay(sub_df,
                         time_col = "Time",
                         expression_col = "expression",
                         method = "nls")

        # Fallback to LM if NLS fails
        if (is.na(res$coef$k)) {
            res <- fit_decay(sub_df,
                             time_col = "Time",
                             expression_col = "expression",
                             method = "lm")
        }

        k_val <- res$coef$k
        half_life <- NA
        if (!is.na(k_val) && k_val > 0) {
            half_life <- log(2) / k_val
        }

        res_rows[[cond]] <- data.frame(
            Geneid = gene_id,
            Condition = cond,
            k = k_val,
            k_se = res$coef$k_se,
            half_life = half_life,
            method = res$coef$method,
            stringsAsFactors = FALSE
        )
    }

    return(do.call(rbind, res_rows))
}

# Run with progress reporting
all_results <- vector("list", n_genes)

for (i in 1:n_genes) {
    if (i %% 500 == 0) print(paste("Processed", i, "/", n_genes, "genes"))
    all_results[[i]] <- process_gene(i)
}

# ---- Combine and annotate results ----
final_results <- do.call(rbind, all_results)

# Merge motif information
motif_info <- data[, c("Geneid", "motif_count_UAAAG", "motif_count_UAAGA",
                       "total_motifs", "has_motif")]
final_results <- merge(final_results, motif_info, by = "Geneid", all.x = TRUE)

# ---- Write output ----
print(paste("Writing results to:", output_file))
write.csv(final_results, output_file, row.names = FALSE)
print("Done.")
