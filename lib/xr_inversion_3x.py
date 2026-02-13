#!/usr/bin/env python3

import pandas as pd
import os

# ============================
# Hardcoded input CSVs
# ============================
MODEL1_CSV = "/home/marchandlab/DataAnalysis/Kaplan/basecall/DsPx/251217_D10_BC/TPxC-N-BC/demux/demux_per-read_modifications.tsv"  # Model 1: standard
MODEL2_CSV = "/home/marchandlab/DataAnalysis/Kaplan/basecall/DsPx/251217_D10_BC/TDsC-N-BC/demux/demux_per-read_modifications.tsv"  # Model 2: RC 
MODEL3_CSV = "/home/marchandlab/DataAnalysis/Kaplan/basecall/DsPx/251217_D10_BC/TDsC-Px-BC/demux/demux_per-read_modifications.tsv"  # Model 3: Px vs Ds




OUTDIR = "/home/marchandlab/DataAnalysis/Kaplan/basecall/DsPx/251217_D10_BC/TPxC-Inversions"
os.makedirs(OUTDIR, exist_ok=True)

PER_READ_OUT = os.path.join(OUTDIR, "per_read_comparison.csv")
SUMMARY_OUT = os.path.join(OUTDIR, "per_barcode_pair_summary.csv")


# ============================
# Load input TSVs
# ============================
m1 = pd.read_csv(MODEL1_CSV, sep="\t")
m2 = pd.read_csv(MODEL2_CSV, sep="\t")
m3 = pd.read_csv(MODEL3_CSV, sep="\t")

# Ensure read_id is string
for df in (m1, m2, m3):
    df["read_id"] = df["read_id"].astype(str)


# ============================
# Merge all three models
# ============================
merged12 = m1.merge(
    m2,
    on="read_id",
    suffixes=("_m1", "_m2"),
    how="inner"
)

merged = merged12.merge(
    m3,
    on="read_id",
    suffixes=("", "_m3"),   # model3 gets "_m3"
    how="inner"
)

# For Model 3 columns, manually append suffix
merged = merged.rename(
    columns={
        "class_pred": "class_pred_m3",
        "class_probs": "class_probs_m3",
        "barcode_pair": "barcode_pair_m3"
    }
)


# ============================
# Build per-read dataframe
# ============================
per_read_df = merged[
    [
        "read_id",
        "barcode_pair_m1",     # using m1 as canonical barcode
        "class_pred_m1",
        "class_probs_m1",
        "class_pred_m2",
        "class_probs_m2",
        "class_pred_m3",
        "class_probs_m3",
    ]
].rename(columns={"barcode_pair_m1": "barcode_pair"})


# ============================
# Save per-read comparison file
# ============================
per_read_df.to_csv(PER_READ_OUT, index=False)


# ============================
# Compute per-barcode summary
# ============================

def compare_3(row):
    """
    Returns a 3-character code:
      e.g. "010", "111", "001", ...
    representing (m1, m2, m3) predictions.
    """
    m1 = int(row["class_pred_m1"])
    m2 = int(row["class_pred_m2"])
    m3 = int(row["class_pred_m3"])
    return f"{m1}{m2}{m3}"

per_read_df["comparison"] = per_read_df.apply(compare_3, axis=1)

# Create indicator columns for all 8 possible patterns
patterns = [f"{a}{b}{c}" for a in "01" for b in "01" for c in "01"]

for p in patterns:
    per_read_df[f"pat_{p}"] = (per_read_df["comparison"] == p).astype(int)

# Group summary by barcode pair
summary = (
    per_read_df
    .groupby("barcode_pair")
    .agg(
        total_reads=("read_id", "count"),
        **{f"frac_{p}": (f"pat_{p}", "mean") for p in patterns}
    )
    .reset_index()
)

summary.to_csv(SUMMARY_OUT, index=False)

print("Finished!")
print(f"Per-read file saved to: {PER_READ_OUT}")
print(f"Summary file saved to: {SUMMARY_OUT}")

