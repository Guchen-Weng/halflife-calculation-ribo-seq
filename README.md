# halflife-calculation-ribo-seq

RNA-seq half-life and Ribo-seq analysis pipelines for *E. coli* K-12 MG1655 (EV vs WT).

## Environment

### Conda (recommended)

```bash
conda create -n EcoliRNA python=3.10
conda activate EcoliRNA
conda install -c bioconda fastp bowtie2 samtools subread bedtools
pip install pandas numpy seaborn matplotlib
```

### R packages

```r
# CRAN
install.packages(c("ggplot2", "dplyr", "tidyr", "ggrepel", "cowplot",
                   "pheatmap", "RColorBrewer", "scales", "data.table",
                   "patchwork", "stringr", "tibble"))

# Bioconductor
BiocManager::install(c("DESeq2", "clusterProfiler", "enrichplot",
                       "org.EcK12.eg.db", "GenomicRanges", "rtracklayer",
                       "Biostrings"))
```

### Tools

| Tool | Purpose |
|------|---------|
| fastp ≥ 0.23 | Read QC |
| Bowtie2 ≥ 2.5 | Alignment |
| samtools ≥ 1.17 | SAM/BAM processing |
| featureCounts | Gene quantification |
| bedtools | BED operations |

## Workflow

### RNA-seq half-life pipeline

Run from `Halflife_analysis/`:

1. `bash 00_fastp.bash` — Raw read QC
2. `bash 01-alignment-double.bash` — Genome alignment
3. `python 02-Readcount.py` — Merge count matrices
4. `Rscript 03_Normalization/tpm_boxplot.R` — TPM normalization
5. `Rscript 04_halflife/calculate_half_life.R` — Half-life estimation
6. `Rscript 04_halflife/plot_half_life_*.R` — Visualization
7. `cd deseq/ && Rscript run_deseq2.R` — Differential expression + GO

### Reference index

Before alignment, build the Bowtie2 genome index from each pipeline's `reference/` directory:

```bash
cd reference/ && bash build_index.sh
### Ribo-seq pipeline

Run from `Riboseq_analysis/`:

1. `bash 00_fastp.bash` — Raw read QC
2. `bash run_filter.sh` — tRNA/rRNA removal
3. `bash 02_alignment.sh` — Genome alignment
4. `bash 03_TPM/run_featurecounts_tpm.sh && python 03_TPM/calculate_tpm.py` — TPM
5. `bedtools bamtobed -i <sample>.sorted.bam > <sample>.bed` — BAM → BED
6. `Rscript 04_monosome/01_annotate_reads.R` — CDS annotation
7. `Rscript 04_monosome/02_pause_score_calculation.R` — Pause scores
8. `Rscript 04_monosome/03_plot_pause_scores.R` — EV vs WT plots
9. `Rscript 05_disome/build_disome_cache_from_endpoint_tables.R` — Disome cache
10. `Rscript 05_disome/plot_disome_schemeA_heatmaps.R` — Disome heatmaps


```

