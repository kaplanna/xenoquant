#!/usr/bin/env python3

import pandas as pd
import os

# ============================
# Hardcoded input CSVs
# ============================
MODEL1_CSV = "/home/marchandlab/DataAnalysis/Kaplan/basecall/8letter/251013_P8_H4_Basecall/SN-Basecall/demux/demux_per-read_modifications.tsv" #Model 1 BN / SN
MODEL2_CSV = "/home/marchandlab/DataAnalysis/Kaplan/basecall/8letter/251013_P8_H4_Basecall/ZC-Basecall/demux/demux_per-read_modifications.tsv" #Model 2 PG / ZN

OUTDIR = "/home/marchandlab/DataAnalysis/Kaplan/basecall/8letter/251013_P8_H4_Basecall/H4-SZ-Combined_Calc"
os.makedirs(OUTDIR, exist_ok=True)

PER_READ_OUT = os.path.join(OUTDIR, "per_read_comparison.csv")
SUMMARY_OUT = os.path.join(OUTDIR, "per_barcode_pair_summary.csv")


# ============================
# Load CSVs
# ============================
m1 = pd.read_csv(MODEL1_CSV, sep="\t")
m2 = pd.read_csv(MODEL2_CSV, sep="\t")

# Make sure read_id is treated as string
m1["read_id"] = m1["read_id"].astype(str)
m2["read_id"] = m2["read_id"].astype(str)


# ============================
# Merge on read_id
# ============================
merged = m1.merge(
    m2,
    on="read_id",
    suffixes=("_m1", "_m2"),
    how="inner"  # only reads found in both files
)

# Keep relevant columns
per_read_df = merged[
    [
        "read_id",
        "barcode_pair_m1",        # same as model2 ideally but using m1
        "class_pred_m1",
        "class_probs_m1",
        "class_pred_m2",
        "class_probs_m2",
    ]
]

# Rename barcode pair for clarity
per_read_df = per_read_df.rename(columns={"barcode_pair_m1": "barcode_pair"})


# ============================
# Save per-read comparison file
# ============================
per_read_df.to_csv(PER_READ_OUT, index=False)


# ============================
# Compute per-barcode summary
# ============================

def compare_classes(row):
    """
    Returns comparison category:
      - "00"   → both models predict 0
      - "11"   → both predict 1
      - "01"   → model1=0, model2=1
      - "10"   → model1=1, model2=0
    """
    m1 = int(row["class_pred_m1"])
    m2 = int(row["class_pred_m2"])
    return f"{m1}{m2}"

per_read_df["comparison"] = merged.apply(compare_classes, axis=1)

# For readability, create label columns
per_read_df["match_00"] = (per_read_df["comparison"] == "00").astype(int)
per_read_df["match_11"] = (per_read_df["comparison"] == "11").astype(int)
per_read_df["diff_01"] = (per_read_df["comparison"] == "01").astype(int)
per_read_df["diff_10"] = (per_read_df["comparison"] == "10").astype(int)


# Group by barcode pair and compute fractions
summary = (
    per_read_df
    .groupby("barcode_pair")
    .agg(
        total_reads=("read_id", "count"),
        frac_00=("match_00", "mean"),
        frac_11=("match_11", "mean"),
        frac_01=("diff_01", "mean"),
        frac_10=("diff_10", "mean"),
    )
    .reset_index()
)

summary.to_csv(SUMMARY_OUT, index=False)

print("Finished!")
print(f"Per-read file saved to: {PER_READ_OUT}")
print(f"Summary file saved to: {SUMMARY_OUT}")

