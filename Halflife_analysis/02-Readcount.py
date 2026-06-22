# ============================================================
# Script:    02-Readcount.py
# Pipeline:  Step 02 - Gene Count Merging
# Purpose:   Merge mRNA and ncRNA gene count matrices by finding
#            common columns and concatenating vertically. Outputs
#            a single combined count matrix for downstream analysis.
# Input:     ../02-Readcount/merged_gene_counts_with_annotations.csv
#            ../02-Readcount/merged_gene_counts_with_annotations_ncRNA.csv
# Output:    ../02-Readcount/merged_gene_counts_final.csv
# Author:    Guchen-Weng
# Date:      2025-01
# ============================================================

import pandas as pd
import sys

# ---- Load data ----
df1 = pd.read_csv("../02-Readcount/merged_gene_counts_with_annotations.csv")
df2 = pd.read_csv("../02-Readcount/merged_gene_counts_with_annotations_ncRNA.csv")

# ---- Merge ----
# Find common columns (metadata + shared sample columns)
common_cols = list(set(df1.columns) & set(df2.columns))

# Vertical concatenation using shared columns
df_merged = pd.concat([df1[common_cols], df2[common_cols]], ignore_index=True)

# ---- Restore column order: metadata first, then samples ----
base_cols = ['Geneid', 'chromosome', 'start', 'end', 'strand',
             'biotype', 'description', 'gene_name']
sample_cols = sorted([c for c in df_merged.columns if c not in base_cols])
final_cols = base_cols + sample_cols

df_final = df_merged[final_cols]

# ---- Write output ----
output_file = "../02-Readcount/merged_gene_counts_final.csv"
df_final.to_csv(output_file, index=False)
print(f"Merged file saved to {output_file}")
