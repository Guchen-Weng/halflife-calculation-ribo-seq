# Reference Data

All genome files obtained from Ensembl Genomes (release 62) for *Escherichia coli* str. K-12 substr. MG1655 (assembly ASM584v2).

## Genome & Annotation

| File | Description | Size |
|------|-------------|------|
| `Escherichia_coli_...dna.chromosome.Chromosome.fa` | Chromosomal DNA sequence | 4.6 MB |
| `Escherichia_coli_...Chromosome.gff3` | Gene annotation (GFF3) | 2.6 MB |
| `Escherichia_coli_...Chromosome.gtf` | Gene annotation (GTF) | 1.3 MB |
| `Escherichia_coli_...cds.all.fa` | CDS sequences | 4.7 MB |

## Filtering References

| File | Description | Size |
|------|-------------|------|
| `tRNA.fa` | E. coli tRNA sequences | 19 KB |
| `rRNA.fa` | E. coli rRNA sequences | 35 KB |

## Building Indices

### Bowtie2 genome index
```bash
cd reference/
bash build_index.sh
```

### tRNA/rRNA indices
Built automatically by `run_filter.sh` on first run:
```bash
bowtie2-build reference/tRNA.fa reference/tRNA
bowtie2-build reference/rRNA.fa reference/rRNA
```

Requires [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/) >= 2.5.

Index files (`*.bt2`) are excluded from version control by `.gitignore`.
