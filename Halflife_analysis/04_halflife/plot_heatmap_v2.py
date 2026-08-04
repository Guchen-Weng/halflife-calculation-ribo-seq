# ============================================================
# Script:    plot_heatmap_v2.py
# Pipeline:  Step 04 - Half-life Visualization
# Purpose:   Heatmap of residual RNA (TPM normalized to t=0min)
#            for the top 50 genes with largest WT vs. EV half-life
#            difference. Outputs separate heatmaps for WT and EV.
# Input:     ../../02-Readcount/merged_gene_counts_final.csv
#            ../../half_decay/half_life_results.csv
# Output:    ../../03-Visualization/WT_top50_halflife_diff_heatmap.pdf
#            ../../03-Visualization/EV_top50_halflife_diff_heatmap.pdf
# Author:    Guchen-Weng
# ============================================================

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
import re

# Set style
sns.set_theme(style="whitegrid")

# ---- File paths ----
input_file = '../../02-Readcount/merged_gene_counts_final.csv'
half_life_file = '../../half_decay/half_life_results.csv'
output_dir = '../../03-Visualization'
os.makedirs(output_dir, exist_ok=True)

# ---- 1. Load and process half-life data ----
print(f"Loading half-life data from {half_life_file}...")
hl_df = pd.read_csv(half_life_file)

# Pivot to get WT and EV half-lives side by side
hl_pivot = hl_df.pivot(index='Geneid', columns='Condition', values='half_life')

# Keep genes with both WT and EV estimates
hl_filtered = hl_pivot.dropna(subset=['WT', 'EV'])
print(f"Genes with both WT and EV half-lives: {len(hl_filtered)}")

# Select top 50 genes by absolute half-life difference
hl_filtered['diff'] = (hl_filtered['WT'] - hl_filtered['EV']).abs()
top50_genes_df = hl_filtered.sort_values('diff', ascending=False).head(50)
top50_gene_ids = top50_genes_df.index.tolist()

print(f"Selected {len(top50_gene_ids)} genes based on half-life difference.")
print("Top 5 genes:", top50_gene_ids[:5])

# ---- 2. Load expression data ----
print(f"Loading expression data from {input_file}...")
df = pd.read_csv(input_file)

# Calculate gene length
df['start'] = pd.to_numeric(df['start'], errors='coerce')
df['end'] = pd.to_numeric(df['end'], errors='coerce')
df['length'] = (df['end'] - df['start']).abs()

# Identify metadata vs. sample columns
metadata_cols = ['Geneid', 'chromosome', 'start', 'end', 'strand',
                 'biotype', 'description', 'gene_name', 'length']
sample_cols = [c for c in df.columns if c not in metadata_cols]

# ---- 3. Calculate TPM ----
rpk = df[sample_cols].div(df['length'] / 1000, axis=0)
scaling_factor = rpk.sum(axis=0) / 1e6
tpm = rpk.div(scaling_factor, axis=1)

tpm['Geneid'] = df['Geneid']
tpm['gene_name'] = df['gene_name']
tpm = tpm.set_index('Geneid')

# Verify top 50 genes present in TPM data
missing_genes = [g for g in top50_gene_ids if g not in tpm.index]
if missing_genes:
    print(f"Warning: {len(missing_genes)} genes from top 50 "
          f"not found in expression data.")
    print(missing_genes)
else:
    print("All top 50 genes found in expression data.")

# ---- 4. Parse sample info and compute group means ----
sample_info = []
for col in sample_cols:
    new_col = col.replace('omin', '0min')
    parts = new_col.split('_')
    if len(parts) >= 3:
        group = parts[1]       # WT or EV
        time_str = parts[2]    # 0min, 2min, 6min, 12min
        match = re.search(r'(\d+)min', time_str)
        if match:
            time_val = int(match.group(1))
            sample_info.append({
                'col': col,
                'group': group,
                'time': time_val
            })

sample_info_df = pd.DataFrame(sample_info)

# Mean TPM per Group x Time combination
grouped_tpm = pd.DataFrame(index=tpm.index)
grouped_tpm['gene_name'] = tpm['gene_name']

groups = ['WT', 'EV']
times = [0, 2, 6, 12]

for group in groups:
    for time in times:
        samples = sample_info_df[
            (sample_info_df['group'] == group) &
            (sample_info_df['time'] == time)
        ]['col'].tolist()

        if samples:
            col_name = f"{group}_{time}min"
            grouped_tpm[col_name] = tpm[samples].mean(axis=1)
        else:
            print(f"Warning: No samples found for {group} {time}min")

# ---- 5. Normalize to t=0min (residual RNA) ----
normalized_tpm = pd.DataFrame(index=grouped_tpm.index)
normalized_tpm['gene_name'] = grouped_tpm['gene_name']

for group in groups:
    base_col = f"{group}_0min"
    if base_col not in grouped_tpm.columns:
        continue

    base_val = grouped_tpm[base_col]

    # Handle zero expression at t=0 (replace with NaN to force NaN ratio)
    zero_base_count = (base_val == 0).sum()
    if zero_base_count > 0:
        print(f"Warning: {zero_base_count} genes have 0 expression "
              f"at {base_col}.")
        base_val = base_val.replace(0, np.nan)

    for time in times:
        col = f"{group}_{time}min"
        if col in grouped_tpm.columns:
            rel_exp = grouped_tpm[col] / base_val
            normalized_tpm[col] = rel_exp

# ---- 6. Filter for top 50 genes ----
plot_data = normalized_tpm[
    normalized_tpm.index.isin(top50_gene_ids)
].copy()

# Add display name (gene_name or Geneid if missing)
plot_data['display_name'] = plot_data.apply(
    lambda x: x['gene_name'] if pd.notna(x['gene_name']) else x.name, axis=1
)
plot_data = plot_data.reindex(top50_gene_ids)

# Check for NaN values
print("Checking for NaNs in plot data...")
nan_counts = plot_data.isna().sum()
print(nan_counts)

# Fill remaining NaN with 0 for heatmap display
plot_data_filled = plot_data.fillna(0)

# ---- 7. Plot heatmaps ----
def plot_heatmap_v2(group_name, output_filename):
    """Generate heatmap for a single condition group."""
    cols = [f"{group_name}_{t}min" for t in times]
    heatmap_data = plot_data_filled[cols]

    # Rename columns for cleaner display
    heatmap_data.columns = [f"{t}min" for t in times]
    heatmap_data.index = plot_data_filled['display_name']

    plt.figure(figsize=(6, 12))

    ax = sns.heatmap(
        heatmap_data,
        cmap="YlGnBu",
        xticklabels=True,
        yticklabels=True,
        vmin=0,
        vmax=1,
        cbar_kws={'label': 'Residual RNA'}
    )

    plt.title(f'{group_name} Top 50 Genes (Half-life Diff)', fontsize=14)
    plt.xlabel('Time', fontsize=12)
    plt.ylabel('Genes', fontsize=12)
    plt.tight_layout()

    save_path = os.path.join(output_dir, output_filename)
    plt.savefig(save_path, format='pdf', bbox_inches='tight')
    plt.close()
    print(f"Saved heatmap to {save_path}")

# Generate WT and EV heatmaps
plot_heatmap_v2('WT', 'WT_top50_halflife_diff_heatmap.pdf')
plot_heatmap_v2('EV', 'EV_top50_halflife_diff_heatmap.pdf')
