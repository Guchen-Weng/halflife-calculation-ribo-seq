suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

# ===== 1) Working directory and input files =====
workdir <- "/Users/fanyanglv/Riboseq0127"
if (!dir.exists(workdir)) stop("Directory not found: ", workdir)
setwd(workdir)
message("Working directory: ", getwd())

files <- c(
  "EV_1_cds_annotation.csv",
  "EV_2_cds_annotation.csv",
  "WT_1_cds_annotation.csv",
  "WT_2_cds_annotation.csv"
)
stopifnot(all(file.exists(files)))

# ===== 2) Read-length filtering range =====
len_min <- 26
len_max <- 32

# ===== 3) Batch filtering and read-count summary =====
stats <- rbindlist(lapply(files, function(f) {
  dt <- fread(f)
  
  # Convert relevant columns to numeric types without applying additional row filtering
  dt[, read_length := as.integer(read_length)]
  dt[, readcount   := as.numeric(readcount)]
  
  # Total read counts before read-length filtering
  reads_before <- dt[, sum(readcount, na.rm = TRUE)]
  
  # Retain reads with lengths between len_min and len_max, inclusive
  dt_filt <- dt[read_length >= len_min & read_length <= len_max]
  reads_after <- dt_filt[, sum(readcount, na.rm = TRUE)]
  
  # Write the filtered table to a new CSV file
  out_file <- str_replace(
    f,
    "\\.csv$",
    paste0("_len", len_min, "-", len_max, ".csv")
  )
  fwrite(dt_filt, out_file)
  
  # Return per-sample filtering statistics
  data.table(
    sample = str_replace(f, "_cds_annotation\\.csv$", ""),
    in_file = f,
    out_file = out_file,
    reads_before = reads_before,
    reads_after  = reads_after,
    frac_kept = ifelse(reads_before > 0, reads_after / reads_before, NA_real_)
  )
}), fill = TRUE)

# Save the summary table across all samples
summary_file <- paste0("reads_summary_len", len_min, "-", len_max, ".csv")
fwrite(stats, summary_file)

print(stats)

message("Done. Filtered CSV files and summary table were saved in: ", getwd())
message("Summary file: ", summary_file)