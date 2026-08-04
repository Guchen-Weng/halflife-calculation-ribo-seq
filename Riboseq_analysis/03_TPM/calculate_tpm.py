# ============================================================
# Script:    calculate_tpm.py
# Pipeline:  Step 03 - TPM Calculation
# Purpose:   Calculate TPM (Transcripts Per Million) from
#            featureCounts output
# Input:     03_TPM/counts.txt (featureCounts output)
# Output:    03_TPM/tpm_counts.csv
# Author:    Guchen-Weng
# ============================================================

import pandas as pd
import sys
import os

if not os.path.exists("counts.txt"):
    print("Error: counts.txt not found!")
    sys.exit(1)

# featureCounts output has a comment line (cmd) starting with #
try:
    df = pd.read_csv("counts.txt", sep="\t", comment="#")
except Exception as e:
    print(f"Error reading counts.txt: {e}")
    sys.exit(1)

# Column layout: Geneid, Chr, Start, End, Strand, Length, Sample1, ...
sample_cols = df.columns[6:]

print(f"Calculating TPM for samples: {list(sample_cols)}")

# RPK = reads per kilobase
lengths_kb = df["Length"] / 1000
rpk = df[sample_cols].div(lengths_kb, axis=0)

# TPM = RPK / (sum(RPK) / 1e6)
scaling_factor = rpk.sum(axis=0) / 1e6
tpm = rpk.div(scaling_factor, axis=1)

tpm.insert(0, "Geneid", df["Geneid"])

output_file = "tpm_counts.csv"
tpm.to_csv(output_file, index=False)
print(f"TPM calculation complete. Saved to {output_file}")
