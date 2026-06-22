# Ribo-seq Analysis Pipeline (Monosome + Disome)

Analysis pipeline for Ribo-seq data in *E. coli* K-12, comparing EV (empty vector) and WT (wild-type) conditions. Covers raw read QC → alignment → TPM → monosome pause scores → disome endpoint analysis.

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
├── 04_monosome/
│   ├── 01_annotate_reads.R               # Annotate reads to CDS via GFF3
│   ├── 02_pause_score_calculation.R      # Per-codon/AA/asymmetry scores
│   └── 03_plot_pause_scores.R            # EV vs WT pause score plots
├── 05_disome/
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
       ├──► bedtools bamtobed    BAM → BED conversion
       │         │
       │    04_monosome/
       │    ├── 01_annotate_reads.R       GFF3 → CDS annotation of reads
       │    ├── 02_pause_score_calculation.R   Codon/AA/Asymmetry scores
       │    └── 03_plot_pause_scores.R         EV vs WT scatter + log2FC plots
       │
       └──► 05_disome/           (Read length 35-65 nt)
            ├── build_disome_cache...R   z-scores for stop codons + 5 motifs
            └── plot_disome_schemeA_heatmaps.R  dZ (WT-EV) heatmaps + mean curves
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
| [bedtools](https://bedtools.readthedocs.io/) | BAM to BED conversion + BED operations |

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
| 4a | `bedtools bamtobed -i <sample>.sorted.bam > <sample>.bed` | `../02-alignment/` | `04_monosome/genome-aligned-bam/` |
| 4b | `Rscript 04_monosome/01_annotate_reads.R` | reference/ GFF + BED files | `04_monosome/genome-aligned-bam/` |
| 4c | `Rscript 04_monosome/02_pause_score_calculation.R` | `*_with_PEA_offset15_fixedStrand.csv` + reference/ cds.all.fa | `04_monosome/PauseScore_fixedStrand_len26-32/` |
| 4d | `Rscript 04_monosome/03_plot_pause_scores.R` | PauseScore CSV files | EV vs WT plots |
| 5a | `Rscript 05_disome/build_disome_cache_from_endpoint_tables.R` | genome-aligned-bam/ + reference/ | `05_disome/cache_stop_and_motifs_hm_dz.rds` |
| 5b | `Rscript 05_disome/plot_disome_schemeA_heatmaps.R` | cache .rds | `05_disome/SchemeA_*/` |

### BAM to BED conversion (Step 4a)

After genome alignment, convert sorted BAM files to BED format before annotation. Run from `Riboseq_analysis/` root:

```bash
# Create output directory
mkdir -p 04_monosome/genome-aligned-bam

# Convert each sample
for bam in ../02-alignment/*.sorted.bam; do
    name=$(basename "$bam" .sorted.bam)
    bedtools bamtobed -i "$bam" > "04_monosome/genome-aligned-bam/${name}.bed"
done
```

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
- **Pause score calculation** (Step 4c): Uses pre-computed P/A/E-site positions with offset 15. Input files (`*_with_PEA_offset15_fixedStrand.csv`) are produced by a prior annotation step that adds A/P/E-site coordinates relative to CDS start.
- **Shared upstream**: Steps 0-3 produce data for both this pipeline and [Halflife_analysis](../Halflife_analysis/). Run steps 0-3 once, then proceed with either or both downstream pipelines.
