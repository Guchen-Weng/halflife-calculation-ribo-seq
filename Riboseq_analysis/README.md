# Ribo-seq Analysis Pipeline (Monosome + Disome)

Analysis pipeline for Ribo-seq data in *E. coli* K-12, comparing EV (empty vector) and WT (wild-type) conditions. Covers raw read QC → alignment → TPM → CDS processing → monosome pause scores → disome endpoint analysis.

This pipeline shares upstream steps (QC, alignment, TPM) with the companion [Halflife_analysis](../Halflife_analysis/) pipeline, which performs mRNA half-life estimation from RNA-seq data.

## Directory Structure

```
Riboseq_analysis/
├── 00_fastp.bash                         # Raw FASTQ QC (fastp)
├── run_filter.sh                         # tRNA/rRNA removal (Bowtie2)
├── 02_alignment.sh                       # Bowtie2 genome alignment
├── 03_TPM/
│   ├── run_featurecounts_tpm.sh          # featureCounts gene quantification
│   └── calculate_tpm.py                  # TPM normalization
├── 04_CDS/
│   ├── buildnonoveralpcds.bash           # Build non-overlapping CDS index
│   └── run_cds_align.sh                  # Align reads to CDS
├── 05_monosome/
│   ├── annotate_reads.R                  # Annotate reads to CDS via GFF3
│   ├── calculate_offset_table.R          # P-site offset optimization
│   ├── process_bam_csv.R                 # Filter reads + compute P-site
│   ├── process_trim_overlap.R            # Trim CDS ends + replicate QC
│   ├── process_pause_score_nofilter.R    # Compute codon pause scores
│   ├── 01_pause_score_calculation.R      # Per-codon/AA/asymmetry scores
│   └── 02_plot_pause_scores.R            # EV vs WT pause score plots
├── 06_disome/
│   ├── build_disome_cache_from_endpoint_tables.R  # Disome z-score cache
│   └── plot_disome_schemeA_heatmaps.R             # dZ heatmaps + curves
├── reference/                            # Reference genome + annotations
├── .gitignore
└── README.md
```

## Pipeline Overview

```
00_fastp.bash                    Raw FASTQ QC (fastp)
       │
run_filter.sh                    tRNA/rRNA removal
       │
02_alignment.sh                  Bowtie2 genome alignment
       │
       ├──► 03_TPM/              featureCounts → TPM
       │         │
       │         └──► Halflife_analysis/  (mRNA half-life pipeline)
       │
       └──► 04_CDS/              Build CDS index + align to CDS
                 │
                 ▼
05_monosome/
├── annotate_reads.R             GFF3 → CDS annotation of reads
├── calculate_offset_table.R     Find optimal P-site offset
├── process_bam_csv.R            Filter (len 26-30) + P-site (offset 13)
├── process_trim_overlap.R       Trim CDS ends + replicate overlap check
├── process_pause_score_nofilter.R   Codon pause scores (no zero-fill)
├── 01_pause_score_calculation.R     Codon/AA/Asymmetry (offset 15, len 26-32)
└── 02_plot_pause_scores.R           EV vs WT scatter + log2FC plots
                 │
06_disome/                       (Read length 35-65 nt)
├── build_disome_cache...R       z-scores for stop codons + 5 motifs
└── plot_disome_schemeA_heatmaps.R   dZ (WT-EV) heatmaps + mean curves
```

## Requirements

Shared with [Halflife_analysis](../Halflife_analysis/). See its README for full dependency list.

### R packages (additional)

```r
# CRAN
install.packages(c("data.table", "patchwork", "scales"))

# Bioconductor
BiocManager::install(c("GenomicRanges", "rtracklayer", "Biostrings"))
```

### Command-line tools (additional)

| Tool | Purpose |
|------|---------|
| [featureCounts](http://subread.sourceforge.net/) | Gene-level read counting |
| [bedtools](https://bedtools.readthedocs.io/) | BED file operations |

### Conda environment (shared)

```bash
conda create -n EcoliRNA python=3.10
conda activate EcoliRNA
conda install -c bioconda fastp bowtie2 samtools subread bedtools
pip install pandas numpy seaborn matplotlib
```

## Execution Order

All scripts run from `Riboseq_analysis/` root.

| Step | Script | Input | Output |
|------|--------|-------|--------|
| 0 | `bash 00_fastp.bash` | `../RawData/` | `../00-QC/` |
| 1 | `bash run_filter.sh` | `../00-QC/cleaned/` | `../01-filter/` |
| 2 | `bash 02_alignment.sh` | `../01-filter/` | `../02-alignment/` |
| 3 | `bash 03_TPM/run_featurecounts_tpm.sh` | `../02-alignment/` | `03_TPM/` |
| 4a | `bash 04_CDS/buildnonoveralpcds.bash` | reference/ GFF + FASTA | `../03.1-build-non-overlapping-cds/` |
| 4b | `bash 04_CDS/run_cds_align.sh` | `../01-filter/` + CDS index | `../03.2-alignmenttocds/` |
| 5a | `Rscript 05_monosome/annotate_reads.R` | reference/ GFF + BED files | `05_monosome/genome-aligned-bam/` |
| 5b | `Rscript 05_monosome/calculate_offset_table.R` | genome-aligned-bam/ | `05_monosome/offset_optimization/` |
| 5c | `Rscript 05_monosome/process_bam_csv.R` | genome-aligned-bam/ | `05_monosome/filter-26-30-off13/` |
| 5d | `Rscript 05_monosome/process_trim_overlap.R` | filter-26-30-off13/ | `05_monosome/filter-26-30-off13/trim/` |
| 5e | `Rscript 05_monosome/process_pause_score_nofilter.R` | filter + trim + reference/ FASTA | `05_monosome/filter-26-30-off13/pausescore-0nofilter/` |
| 5f* | `Rscript 05_monosome/01_pause_score_calculation.R` | PEA_offset15 CSV + reference/ cds.all.fa | `05_monosome/PauseScore_fixedStrand_len26-32/` |
| 5g* | `Rscript 05_monosome/02_plot_pause_scores.R` | PauseScore CSV files | EV vs WT plots |
| 6a | `Rscript 06_disome/build_disome_cache...R` | genome-aligned-bam/ + reference/ | `06_disome/cache_stop_and_motifs_hm_dz.rds` |
| 6b | `Rscript 06_disome/plot_disome_schemeA_heatmaps.R` | cache .rds | `06_disome/SchemeA_*/` |

*Steps 5f-5g require `*_cds_annotation_with_PEA_offset15_fixedStrand.csv` input files (see Notes).

## Reference Data

All reference files are provided in `reference/`. See `reference/README.md` for sources.

Before running Step 2, build the Bowtie2 genome index:
```bash
cd reference/ && bash build_index.sh
```

tRNA/rRNA indices are built automatically by `run_filter.sh`.

## Output Files

All generated figures (PDF/PNG), data files (CSV), and intermediate outputs are excluded from version control by `.gitignore`.

## Notes

- **Sample naming**: `WT_1`, `WT_2`, `EV_1`, `EV_2` (2 biological replicates per condition).
- **Reference genome**: *Escherichia coli* str. K-12 substr. MG1655 (Ensembl release 62).
- **Monosome vs. Disome**: Monosome reads = 26-32 nt (standard Ribo-seq). Disome reads = 35-65 nt (paired ribosomes, distinct biological signal).
- **Pause score calculation** (steps 5f-5g): Uses pre-computed P/A/E-site positions with offset 15. The input files (`*_with_PEA_offset15_fixedStrand.csv`) are produced by a prior annotation step that adds A/P/E-site coordinates relative to CDS start. The offset table from step 5b (offset 13 for filter-26-30) uses a different parameterization optimized for the no-filter pause score pipeline.
- **Shared upstream**: Steps 0-3 produce data for both this pipeline and [Halflife_analysis](../Halflife_analysis/). Run steps 0-3 once, then proceed with either or both downstream pipelines.
