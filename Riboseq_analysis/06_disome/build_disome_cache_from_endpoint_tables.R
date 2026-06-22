#!/usr/bin/env Rscript
# ============================================================
# Script:    build_disome_cache_from_endpoint_tables.R
# Pipeline:  Disome Analysis - Cache Building
# Purpose:   Build an RDS cache from disome-seq endpoint tables
#            containing raw endpoint counts, per-sample z-scores,
#            group mean z-scores, and WT-minus-EV dz for stop
#            codons and selected motifs (ATA, TAA, TAAG, CTC, GCT).
#            Read length range: 35-65 nt (disome).
#
#            z is computed independently for every sample, read
#            length, and endpoint window. Therefore dz = WT - EV
#            describes a DIFFERENCE IN LOCAL SHAPE of endpoint
#            distribution, not an absolute coverage difference.
# Input:     05_monosome/genome-aligned-bam/*_cds_annotation.csv
#            reference/Escherichia_coli_...Chromosome.gff3
#            reference/Escherichia_coli_...dna.chromosome.Chromosome.fa
# Output:    cache_stop_and_motifs_hm_dz.rds
# Author:    Guchen-Weng
# Date:      2025-01
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
  library(Biostrings)
})

# =====================================================================
# 0) Paths and sample configuration
# =====================================================================
# All paths relative to Riboseq_analysis/ root
gff_file <- "reference/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.62.chromosome.Chromosome.gff3"
genome_fasta <- "reference/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.dna.chromosome.Chromosome.fa"

# CDS-annotated read tables from monosome step
sample_sheet <- data.table(
  sample = c("EV_1", "EV_2", "WT_1", "WT_2"),
  group  = c("EV",   "EV",   "WT",   "WT"),
  file   = file.path("05_monosome", "genome-aligned-bam", c(
    "EV_1_cds_annotation.csv",
    "EV_2_cds_annotation.csv",
    "WT_1_cds_annotation.csv",
    "WT_2_cds_annotation.csv"
  ))
)

# Disome read length range and endpoint windows
read_lengths <- 35:65
d5 <- -60:-10   # 5' endpoint positions relative to the anchor
d3 <- -10:45    # 3' endpoint positions relative to the anchor

out_rds <- "cache_stop_and_motifs_hm_dz.rds"

# Motif settings (anchor_pos is 1-based position within the motif)
motif_cfg <- data.table(
  target     = c("ATA", "TAA", "TAAG", "CTC", "GCT"),
  motif      = c("ATA", "TAA", "TAAG", "CTC", "GCT"),
  anchor_pos = c(3L,    1L,    1L,     2L,    2L),
  is_Acenter = c(TRUE,  TRUE,  TRUE,   FALSE, FALSE)
)

# =====================================================================
# 1) Read and prepare CDS-annotated read tables
# =====================================================================

read_annotated_reads <- function(sample_sheet, read_lengths) {
  required <- c("cds_id", "strand", "cds_5",
                "read_5", "read_3", "read_length", "readcount")

  ans <- rbindlist(lapply(seq_len(nrow(sample_sheet)), function(i) {
    f <- sample_sheet$file[i]
    if (!file.exists(f)) stop("Missing annotated read file: ", f)

    x <- fread(f)
    missing_cols <- setdiff(required, names(x))
    if (length(missing_cols) > 0) {
      stop("Missing required columns in ", f, ": ",
           paste(missing_cols, collapse = ", "))
    }

    x[, sample := sample_sheet$sample[i]]
    x[, `:=`(
      cds_id      = as.character(cds_id),
      strand      = as.character(strand),
      cds_5       = as.integer(cds_5),
      read_5      = as.integer(read_5),
      read_3      = as.integer(read_3),
      read_length = as.integer(read_length),
      readcount   = as.numeric(readcount)
    )]

    # Convert genomic coordinates to transcript-direction coordinates.
    # Values can be <1 or >CDS length for reads crossing a CDS boundary,
    # which is desired for endpoint profiling around the stop codon.
    x[, read5_tx := ifelse(strand == "+",
                           read_5 - cds_5 + 1L,
                           cds_5 - read_5 + 1L)]
    x[, read3_tx := ifelse(strand == "+",
                           read_3 - cds_5 + 1L,
                           cds_5 - read_3 + 1L)]

    x[read_length %in% read_lengths & !is.na(readcount) & readcount > 0,
      .(sample, cds_id, read_length, read5_tx, read3_tx, readcount)]
  }), use.names = TRUE, fill = TRUE)

  ans[]
}

# =====================================================================
# 2) Build transcript-oriented CDS sequences from GFF3 + genome FASTA
# =====================================================================

get_cds_id <- function(cds_gr) {
  if ("ID" %in% names(mcols(cds_gr))) {
    return(as.character(mcols(cds_gr)$ID))
  }
  if ("Parent" %in% names(mcols(cds_gr))) {
    return(as.character(mcols(cds_gr)$Parent))
  }
  stop("The CDS GFF entries have neither ID nor Parent attributes.")
}

make_cds_sequences <- function(gff_file, genome_fasta) {
  cds_gr <- import.gff3(gff_file)
  cds_gr <- cds_gr[cds_gr$type == "CDS"]
  cds_id <- get_cds_id(cds_gr)

  if (anyNA(cds_id) || any(cds_id == "")) {
    stop("Some CDS entries do not have valid IDs.")
  }

  if (anyDuplicated(cds_id)) {
    stop("Duplicated CDS IDs were detected. This script assumes one contiguous ",
         "CDS record per ID; concatenate exons per CDS before continuing.")
  }

  genome <- readDNAStringSet(genome_fasta)
  seqs <- vapply(seq_along(cds_gr), function(i) {
    chr <- as.character(seqnames(cds_gr)[i])
    idx <- match(chr, names(genome))
    if (is.na(idx)) {
      stop("GFF chromosome name '", chr,
           "' does not match any FASTA sequence name.")
    }

    s <- subseq(
      genome[[idx]],
      start = start(cds_gr)[i],
      end   = end(cds_gr)[i]
    )
    if (as.character(strand(cds_gr)[i]) == "-") {
      s <- reverseComplement(s)
    }
    as.character(s)
  }, character(1))

  cds_seq <- DNAStringSet(seqs)
  names(cds_seq) <- cds_id

  cds_info <- data.table(
    cds_id     = cds_id,
    cds_length = as.integer(width(cds_gr)),
    strand     = as.character(strand(cds_gr))
  )

  list(cds_info = cds_info, cds_seq = cds_seq)
}

# =====================================================================
# 3) Define stop-codon and motif anchors in transcript coordinates
# =====================================================================

make_stop_anchors <- function(cds_info, cds_seq) {
  # Anchor = first nucleotide of the annotated terminal stop codon.
  # In transcript coordinates, a CDS of length L has stop start at L-2.
  terminal_codon <- vapply(
    seq_along(cds_seq),
    function(i) {
      s <- as.character(cds_seq[[i]])
      if (nchar(s) < 3) return(NA_character_)
      substr(s, nchar(s) - 2L, nchar(s))
    },
    character(1)
  )

  if (any(!terminal_codon %in% c("TAA", "TAG", "TGA"))) {
    warning("Some CDS records do not end in TAA/TAG/TGA. ",
            "Check whether the GFF CDS ranges include stop codons.")
  }

  cds_info[, .(cds_id, anchor_tx = cds_length - 2L)]
}

make_motif_anchors <- function(cds_seq, motif, anchor_pos) {
  hit_index <- vmatchPattern(
    pattern = motif,
    subject = cds_seq,
    fixed   = TRUE
  )
  n_hits <- elementNROWS(hit_index)

  if (!any(n_hits > 0)) {
    return(data.table(cds_id = character(), anchor_tx = integer()))
  }

  data.table(
    cds_id = rep(names(cds_seq), n_hits),
    anchor_tx = as.integer(unlist(start(hit_index))) +
                as.integer(anchor_pos) - 1L
  )
}

# =====================================================================
# 4) Count endpoints around each anchor and calculate z-scores
# =====================================================================

make_anchor_grid <- function(anchors, dists) {
  anchors <- unique(anchors[!is.na(anchor_tx)])

  ans <- anchors[, {
    endpoint_matrix <- outer(anchor_tx, dists, "+")
    data.table(
      endpoint_tx = as.integer(as.vector(endpoint_matrix)),
      dist        = rep(as.integer(dists), each = .N)
    )
  }, by = cds_id]

  ans[]
}

make_heatmap_table <- function(
    reads, anchors, dists, endpoint_column,
    sample_levels, read_lengths
) {
  anchor_grid <- make_anchor_grid(anchors, dists)

  endpoints <- reads[, .(
    count = sum(readcount)
  ), by = .(sample, cds_id, read_length,
            endpoint_tx = get(endpoint_column))]

  if (nrow(anchor_grid) == 0L) {
    observed <- data.table(
      sample = character(), read_length = integer(),
      dist = integer(), count = numeric()
    )
  } else {
    joined <- merge(
      endpoints, anchor_grid,
      by = c("cds_id", "endpoint_tx"),
      all = FALSE, allow.cartesian = TRUE
    )

    observed <- joined[, .(count = sum(count)),
                       by = .(sample, read_length, dist)]
  }

  full_grid <- CJ(
    sample      = as.character(sample_levels),
    read_length = as.integer(read_lengths),
    dist        = as.integer(dists),
    unique = TRUE
  )

  hm <- merge(full_grid, observed,
              by = c("sample", "read_length", "dist"),
              all.x = TRUE, sort = FALSE)
  hm[is.na(count), count := 0]

  # z-score: for each sample and read length, standardize counts
  # over all distances using sample SD (n - 1)
  hm[, z := {
    mu <- mean(count)
    sig <- sd(count)
    if (is.na(sig) || sig == 0) rep(0, .N) else (count - mu) / sig
  }, by = .(sample, read_length)]

  setorder(hm, sample, read_length, dist)
  hm[]
}

make_group_mean_and_dz <- function(hm, sample_sheet) {
  x <- merge(hm, sample_sheet[, .(sample, group)],
             by = "sample", all.x = TRUE, sort = FALSE)

  if (anyNA(x$group)) {
    stop("At least one sample in the heatmap table has no group assignment.")
  }

  mean_dt <- x[, .(z = mean(z)), by = .(group, read_length, dist)]

  wide <- dcast(mean_dt, read_length + dist ~ group, value.var = "z")

  if (!all(c("EV", "WT") %in% names(wide))) {
    stop("Both EV and WT groups are required to calculate dz = WT - EV.")
  }

  dz <- wide[, .(read_length, dist, EV, WT, dz = WT - EV)]

  setorder(mean_dt, group, read_length, dist)
  setorder(dz, read_length, dist)

  list(mean = mean_dt[], dz = dz[])
}

make_target_object <- function(
    reads, anchors, sample_sheet, read_lengths, d5, d3
) {
  hm5 <- make_heatmap_table(
    reads = reads, anchors = anchors, dists = d5,
    endpoint_column = "read5_tx",
    sample_levels = sample_sheet$sample,
    read_lengths = read_lengths
  )

  hm3 <- make_heatmap_table(
    reads = reads, anchors = anchors, dists = d3,
    endpoint_column = "read3_tx",
    sample_levels = sample_sheet$sample,
    read_lengths = read_lengths
  )

  s5 <- make_group_mean_and_dz(hm5, sample_sheet)
  s3 <- make_group_mean_and_dz(hm3, sample_sheet)

  list(hm5 = hm5, hm3 = hm3,
       dz5 = s5$dz, dz3 = s3$dz,
       mean5 = s5$mean, mean3 = s3$mean)
}

# =====================================================================
# 5) Build the cache
# =====================================================================

stopifnot(all(file.exists(sample_sheet$file)))
stopifnot(all(sample_sheet$group %in% c("EV", "WT")))

reads <- read_annotated_reads(sample_sheet, read_lengths)
cds_data <- make_cds_sequences(gff_file, genome_fasta)

n_matched <- sum(reads$cds_id %in% cds_data$cds_info$cds_id)
if (n_matched == 0L) {
  stop("None of the CDS IDs in the read tables match the GFF CDS IDs. ",
       "Use the same GFF ID/Parent field in both steps.")
}
reads <- reads[cds_id %in% cds_data$cds_info$cds_id]

# Stop-codon profile
stop_anchors <- make_stop_anchors(cds_data$cds_info, cds_data$cds_seq)
stop_target <- make_target_object(
  reads, stop_anchors, sample_sheet, read_lengths, d5, d3
)

# Motif profiles
motif_targets <- vector("list", nrow(motif_cfg))
names(motif_targets) <- motif_cfg$target

for (i in seq_len(nrow(motif_cfg))) {
  cfg <- motif_cfg[i]
  message("Building motif target: ", cfg$target)

  anchors <- make_motif_anchors(
    cds_seq = cds_data$cds_seq,
    motif = cfg$motif,
    anchor_pos = cfg$anchor_pos
  )

  target_obj <- make_target_object(
    reads, anchors, sample_sheet, read_lengths, d5, d3
  )

  target_obj$anchor_pos <- as.integer(cfg$anchor_pos)
  target_obj$is_Acenter <- as.logical(cfg$is_Acenter)

  motif_targets[[cfg$target]] <- target_obj
}

out <- list(
  meta = list(
    d5_min = min(d5), d5_max = max(d5),
    d3_min = min(d3), d3_max = max(d3),
    read_lengths = as.integer(read_lengths)
  ),
  stop = stop_target,
  motifs = motif_targets
)

saveRDS(out, out_rds)
message("\nSaved cache to:\n", out_rds)

# =====================================================================
# 6) QA summary
# =====================================================================

message("\nQA summary:")
message("Reads retained: ", nrow(reads))
message("CDS IDs represented in read table: ", uniqueN(reads$cds_id))
message("Stop anchors: ", nrow(stop_anchors))
for (nm in names(motif_targets)) {
  message("Motif anchors [", nm, "]: ",
          nrow(make_motif_anchors(
            cds_data$cds_seq,
            motif_cfg[target == nm, motif],
            motif_cfg[target == nm, anchor_pos]
          )))
}

message("\nExample output structure:")
str(out$stop, max.level = 1)
