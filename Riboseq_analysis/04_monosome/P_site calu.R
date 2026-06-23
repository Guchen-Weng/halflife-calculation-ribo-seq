suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

# =========================
# 0) Working directory and input files
# =========================
workdir <- "/Users/fanyanglv/Ribo0127"   # Update this path as needed
if (!dir.exists(workdir)) stop("Directory not found: ", workdir)
setwd(workdir)
message("WD = ", getwd())

files <- c(
  "EV_1_cds_annotation.csv",
  "EV_2_cds_annotation.csv",
  "WT_1_cds_annotation.csv",
  "WT_2_cds_annotation.csv"
)
missing <- files[!file.exists(files)]
if (length(missing) > 0) {
  stop("Missing file(s):\n", paste(missing, collapse="\n"),
       "\n\nFiles currently available in the working directory:\n",
       paste(list.files(), collapse="\n"))
}

# =========================
# 1) Parameters
# =========================
offset <- 15L  # 3' end to P-site offset
step   <- 3L   # Distance between A/P/E sites: one codon (= 3 nt)

# =========================
# 2) Core function: strand-aware P/A/E site assignment
# =========================
add_PEA_fixed_strand <- function(dt, offset = 15L, step = 3L) {
  
  # ---- Standardize data types ----
  # Note: read_3, cds_5, and cds_3 are genomic coordinates stored as integers
  need_cols <- c("cds_id", "strand", "cds_5", "cds_3", "read_3", "cds_length")
  miss_cols <- setdiff(need_cols, names(dt))
  if (length(miss_cols) > 0) {
    stop("Missing required column(s): ", paste(miss_cols, collapse = ", "))
  }
  
  dt[, `:=`(
    read_3     = as.integer(read_3),
    cds_5      = as.integer(cds_5),
    cds_3      = as.integer(cds_3),
    cds_length = as.integer(cds_length)
  )]
  
  # ---- Define genomic left and right CDS boundaries ----
  # Do not assume cds_5/cds_3 naming is consistently ordered; use coordinate values instead
  dt[, cds_left  := pmin(cds_5, cds_3)]
  dt[, cds_right := pmax(cds_5, cds_3)]
  
  # ---- (A) Calculate genomic P/A/E coordinates in a strand-aware manner ----
  # Plus strand:  P = read_3 - offset
  # Minus strand: P = read_3 + offset
  dt[, P_genome := fifelse(strand == "+", read_3 - offset, read_3 + offset)]
  dt[, A_genome := fifelse(strand == "+", P_genome + step, P_genome - step)]
  dt[, E_genome := fifelse(strand == "+", P_genome - step, P_genome + step)]
  
  # ---- (B) Convert genomic coordinates to CDS coordinates along the translation direction ----
  # Plus strand: coordinates increase from cds_left
  # Minus strand: coordinates decrease from cds_right because translation is opposite to genomic direction
  dt[, P_cds := fifelse(
    strand == "+",
    P_genome - cds_left + 1L,
    cds_right - P_genome + 1L
  )]
  dt[, A_cds := fifelse(
    strand == "+",
    A_genome - cds_left + 1L,
    cds_right - A_genome + 1L
  )]
  dt[, E_cds := fifelse(
    strand == "+",
    E_genome - cds_left + 1L,
    cds_right - E_genome + 1L
  )]
  
  # ---- (C) Determine whether P/A/E sites fall within the CDS ----
  dt[, P_in_cds := (P_cds >= 1L & P_cds <= cds_length)]
  dt[, A_in_cds := (A_cds >= 1L & A_cds <= cds_length)]
  dt[, E_in_cds := (E_cds >= 1L & E_cds <= cds_length)]
  
  # ---- Optional: retain cds_left/cds_right, or remove them if not needed ----
  # dt[, c("cds_left", "cds_right") := NULL]
  
  dt[]
}

# =========================
# 3) Batch processing and output
# =========================
for (f in files) {
  message("Reading: ", f)
  dt <- fread(f)
  
  dt2 <- add_PEA_fixed_strand(dt, offset = offset, step = step)
  
  out <- str_replace(
    f,
    "\\.csv$",
    sprintf("_with_PEA_offset%d_fixedStrand.csv", offset)
  )
  fwrite(dt2, out)
  message("Saved: ", out)
  
  # Quick QC: fraction of P/A/E sites located within CDSs for plus and minus strands
  qc <- dt2[, .(
    n = .N,
    frac_P_in = mean(P_in_cds, na.rm = TRUE),
    frac_A_in = mean(A_in_cds, na.rm = TRUE),
    frac_E_in = mean(E_in_cds, na.rm = TRUE)
  ), by = strand]
  
  message("QC (", f, "):")
  print(qc)
}

message("All done.")