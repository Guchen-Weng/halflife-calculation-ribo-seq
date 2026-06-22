# Reference Data

## Source

All files were obtained from Ensembl Genomes (release 62) for *Escherichia coli* str. K-12 substr. MG1655 (assembly ASM584v2).

| File | Description | Source |
|------|-------------|--------|
| `Escherichia_coli_...dna.chromosome.Chromosome.fa` | Chromosomal DNA sequence (4.6 MB) | [Ensembl FTP](https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-62/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/dna/) |
| `Escherichia_coli_...Chromosome.gff3` | Gene annotation in GFF3 format (2.6 MB) | [Ensembl FTP](https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-62/gff3/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/) |
| `Escherichia_coli_...Chromosome.gtf` | Gene annotation in GTF format (1.3 MB) | Converted from GFF3 |

## Building the Bowtie2 Index

```bash
cd reference/
bash build_index.sh
```

This requires [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/) ≥ 2.5 to be installed and available on PATH.

The index files (`Ecoli.*.bt2`) are excluded from version control by `.gitignore` because they are generated from the FASTA.
