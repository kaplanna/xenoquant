#!/usr/bin/env python3

import pandas as pd
import os

# ============================
# Hardcoded input CSVs
# ============================
MODEL1_CSV = "/home/marchandlab/DataAnalysis/Kaplan/basecall/DsPx/251130_GDsA_Inv_Basecall_Testing/Can/TPxC_TNC/remora_outputs/per-read_modifications.tsv"  # Model 1: standard
MODEL2_CSV = "/home/marchandlab/DataAnalysis/Kaplan/basecall/DsPx/251130_GDsA_Inv_Basecall_Testing/Can/TDsC_TNC/remora_outputs/per-read_modifications.tsv"  # Model 2: RC 
MODEL3_CSV = "/home/marchandlab/DataAnalysis/Kaplan/basecall/DsPx/251130_GDsA_Inv_Basecall_Testing/Can/TDsC_TPxC/remora_outputs/per-read_modifications.tsv"  # Model 3: Px vs Ds

OUTDIR = "/home/marchandlab/DataAnalysis/Kaplan/basecall/DsPx/251130_GDsA_Inv_Basecall_Testing/Can/TPxC_inversion3x"
os.makedirs(OUTDIR, exist_ok=True)

PER_READ_OUT = os.path.join(OUTDIR, "per_read_comparison.csv")
SUMMARY_OUT = os.path.join(OUTDIR, "per_reference_sequence_summary.csv")


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
    suffixes=("", "_m3"),
    how="inner"
)

# Manually rename Model 3 columns
merged = merged.rename(
    columns={
        "class_pred": "class_pred_m3",
        "class_probs": "class_probs_m3",
        "reference_sequence": "reference_sequence_m3"
    }
)


# ============================
# Build per-read dataframe
# ============================
# Using m1 reference_sequence as canonical
per_read_df = merged[
    [
        "read_id",
        "reference_sequence_m1",   # new grouping key
        "class_pred_m1",
        "class_probs_m1",
        "class_pred_m2",
        "class_probs_m2",
        "class_pred_m3",
        "class_probs_m3",
    ]
].rename(columns={"reference_sequence_m1": "reference_sequence"})

per_read_df.to_csv(PER_READ_OUT, index=False)


# ============================
# Compute per-reference summary
# ============================

def compare_3(row):
    """Return a 3-character pattern m1m2m3."""
    m1 = int(row["class_pred_m1"])
    m2 = int(row["class_pred_m2"])
    m3 = int(row["class_pred_m3"])
    return f"{m1}{m2}{m3}"

per_read_df["comparison"] = per_read_df.apply(compare_3, axis=1)

# Make indicator columns for all 8 patterns
patterns = [f"{a}{b}{c}" for a in "01" for b in "01" for c in "01"]

for p in patterns:
    per_read_df[f"pat_{p}"] = (per_read_df["comparison"] == p).astype(int)

summary = (
    per_read_df
    .groupby("reference_sequence")
    .agg(
        total_reads=("read_id", "count"),
        **{f"frac_{p}": (f"pat_{p}", "mean") for p in patterns}
    )
    .reset_index()
)

summary.to_csv(SUMMARY_OUT, index=False)

print("Finished!")
print(f"Per-read file saved to: {PER_READ_OUT}")
print(f"Summary file (grouped by reference_sequence) saved to: {SUMMARY_OUT}")

