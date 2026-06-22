# ============================================================
# Script:    annotate_reads.R
# Pipeline:  Step 05 - Monosome (Ribosome Profiling)
# Purpose:   Annotate BED reads to CDS regions using GFF3 annotation.
#            For each read overlapping a CDS, output the CDS ID,
#            strand, 5' and 3' coordinates (orientation-aware), read
#            length, distance to CDS ends, and aggregated read counts.
# Input:     reference/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.62.chromosome.Chromosome.gff3
#            BED files from 05_monosome/genome-aligned-bam/*.bed
# Output:    05_monosome/genome-aligned-bam/<sample>_cds_annotation.csv
# Author:    Guchen-Weng
# Date:      2025-01
# ============================================================

# ---- Libraries ----
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(tidyverse)
})

# ---- Paths ----
# Run from Riboseq_analysis/ root directory
gff_file <- "reference/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.62.chromosome.Chromosome.gff3"
bed_dir <- "05_monosome/genome-aligned-bam"
output_dir <- "05_monosome/genome-aligned-bam"

# ---- 1. Read and Process GFF ----
message("Reading GFF file...")
gff <- import.gff3(gff_file)

# Filter for CDS features only
cds_db <- gff[gff$type == "CDS"]

# Ensure CDS entries have an ID. Fall back to Parent attribute if needed.
if (is.null(cds_db$ID)) {
  if (!is.null(cds_db$Parent)) {
    cds_db$ID <- as.character(cds_db$Parent)
  } else {
    stop("GFF CDS entries do not have ID or Parent attributes.")
  }
}

message(paste("Loaded", length(cds_db), "CDS regions."))

# ---- 2. Process BED Files ----
bed_files <- list.files(bed_dir, pattern = "\\.bed$", full.names = TRUE)

for (bed_file in bed_files) {
  sample_name <- tools::file_path_sans_ext(basename(bed_file))
  message(paste("Processing sample:", sample_name))

  # Read BED file
  # BED format: chrom, start(0-based), end(1-based exclusive), name, score, strand
  # Convert 0-based start to 1-based: GRanges start = BED start + 1
  bed_df <- read.table(bed_file, header = FALSE, stringsAsFactors = FALSE)

  if (ncol(bed_df) < 6) {
    warning(paste("File", bed_file, "has fewer than 6 columns. Skipping."))
    next
  }

  # Construct GRanges with 1-based coordinates
  reads_gr <- GRanges(
    seqnames = bed_df$V1,
    ranges = IRanges(start = bed_df$V2 + 1, end = bed_df$V3),
    strand = bed_df$V6,
    read_id = bed_df$V4
  )

  # ---- 3. Find Overlaps ----
  # Require strand-specific overlap
  hits <- findOverlaps(reads_gr, cds_db, ignore.strand = FALSE)

  if (length(hits) == 0) {
    message("  No overlaps found.")
    next
  }

  matched_reads <- reads_gr[queryHits(hits)]
  matched_cds <- cds_db[subjectHits(hits)]

  # Extract CDS IDs (flatten lists if nested)
  cds_ids <- matched_cds$ID
  if (is.list(cds_ids)) {
    cds_ids <- sapply(cds_ids, function(x) paste(x, collapse = ";"))
  }

  # Build data frame with orientation-aware coordinates
  # For '+' strand: 5' = start, 3' = end
  # For '-' strand: 5' = end, 3' = start (reverse complement orientation)
  temp_df <- data.frame(
    cds_id = cds_ids,
    strand = as.character(strand(matched_cds)),
    cds_start = start(matched_cds),
    cds_end = end(matched_cds),
    read_start = start(matched_reads),
    read_end = end(matched_reads),
    stringsAsFactors = FALSE
  )

  temp_df <- temp_df %>%
    mutate(
      cds_5 = ifelse(strand == "+", cds_start, cds_end),
      cds_3 = ifelse(strand == "+", cds_end, cds_start),
      read_5 = ifelse(strand == "+", read_start, read_end),
      read_3 = ifelse(strand == "+", read_end, read_start)
    )

  # ---- 4. Aggregate Reads and Annotate ----
  output_df <- temp_df %>%
    group_by(cds_id, strand, cds_5, cds_3, read_5, read_3) %>%
    summarise(readcount = n(), .groups = "drop") %>%
    mutate(
      read_length = abs(read_5 - read_3) + 1,
      cds_length = abs(cds_5 - cds_3) + 1,
      `read3abs(dis_5)` = abs(read_3 - cds_5),
      `read3abs(dis_3)` = abs(read_3 - cds_3)
    ) %>%
    select(cds_id, strand, cds_5, cds_3, read_5, read_3,
           read_length, cds_length, `read3abs(dis_5)`, `read3abs(dis_3)`,
           readcount)

  # ---- 5. Save Output ----
  out_file <- file.path(output_dir, paste0(sample_name, "_cds_annotation.csv"))
  write.csv(output_df, out_file, row.names = FALSE, quote = FALSE)
  message(paste("  Saved results to:", out_file))
}

message("All samples processed.")
