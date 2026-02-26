#!/usr/bin/env python3

import pandas as pd
import os

# ============================
# Hardcoded input CSVs
# ============================
MODEL1_CSV = "/path/to/demux_per-read_modifications.tsv"
MODEL2_CSV = "/path/to/demux_per-read_modifications.tsv"

OUTDIR = "/Output/dir"
os.makedirs(OUTDIR, exist_ok=True)

PER_READ_OUT = os.path.join(OUTDIR, "per_read_comparison.csv")
SUMMARY_OUT = os.path.join(OUTDIR, "per_barcode_pair_summary.csv")

# ============================
# Load TSVs
# ============================
m1 = pd.read_csv(MODEL1_CSV, sep="\t")
m2 = pd.read_csv(MODEL2_CSV, sep="\t")

# Ensure read_id is string
m1["read_id"] = m1["read_id"].astype(str)
m2["read_id"] = m2["read_id"].astype(str)

# ============================
# Track model presence
# ============================
m1["in_model1"] = 1
m2["in_model2"] = 1

# ============================
# Merge on read_id (outer to track missing reads)
# ============================
merged = m1.merge(
    m2,
    on="read_id",
    suffixes=("_m1", "_m2"),
    how="outer"
)

merged["in_model1"] = merged["in_model1"].fillna(0).astype(int)
merged["in_model2"] = merged["in_model2"].fillna(0).astype(int)

merged["called_by_both_models"] = (
    (merged["in_model1"] == 1) & (merged["in_model2"] == 1)
).astype(int)

# Prefer barcode from model 1, otherwise model 2
merged["barcode_pair"] = merged["barcode_pair_m1"].combine_first(
    merged["barcode_pair_m2"]
)

# ============================
# Per-read output (only reads in both models)
# ============================
per_read_df = merged.loc[
    merged["called_by_both_models"] == 1,
    [
        "read_id",
        "barcode_pair",
        "class_pred_m1",
        "class_probs_m1",
        "read_focus_base_m1",
        "class_pred_m2",
        "class_probs_m2",
        "read_focus_base_m2",
        "called_by_both_models"
    ]
].copy()

per_read_df = per_read_df.rename(columns={
    "read_focus_base_m1": "focus_base_m1",
    "read_focus_base_m2": "focus_base_m2"
})

per_read_df.to_csv(PER_READ_OUT, index=False)

# ============================
# Comparison categories
# ============================
def compare_classes(row):
    m1 = int(row["class_pred_m1"])
    m2 = int(row["class_pred_m2"])
    return f"{m1}{m2}"

per_read_df["comparison"] = per_read_df.apply(compare_classes, axis=1)

per_read_df["match_00"] = (per_read_df["comparison"] == "00").astype(int)
per_read_df["match_11"] = (per_read_df["comparison"] == "11").astype(int)
per_read_df["diff_01"]  = (per_read_df["comparison"] == "01").astype(int)
per_read_df["diff_10"]  = (per_read_df["comparison"] == "10").astype(int)

# ============================
# Summary statistics per barcode
# ============================
summary_core = (
    per_read_df
    .groupby("barcode_pair")
    .agg(
        reads_both_models=("read_id", "count"),
        frac_00=("match_00", "mean"),
        frac_11=("match_11", "mean"),
        frac_01=("diff_01", "mean"),
        frac_10=("diff_10", "mean"),
        median_focus_m1=("focus_base_m1", "median"),
        median_focus_m2=("focus_base_m2", "median"),
    )
    .reset_index()
)

# ============================
# Add model-only read counts
# ============================
model1_only = (
    merged
    .loc[(merged["in_model1"] == 1) & (merged["in_model2"] == 0)]
    .groupby("barcode_pair")
    .size()
    .rename("reads_model1_only")
)

model2_only = (
    merged
    .loc[(merged["in_model1"] == 0) & (merged["in_model2"] == 1)]
    .groupby("barcode_pair")
    .size()
    .rename("reads_model2_only")
)

summary = (
    summary_core
    .merge(model1_only, on="barcode_pair", how="left")
    .merge(model2_only, on="barcode_pair", how="left")
)

summary["reads_model1_only"] = summary["reads_model1_only"].fillna(0).astype(int)
summary["reads_model2_only"] = summary["reads_model2_only"].fillna(0).astype(int)

summary.to_csv(SUMMARY_OUT, index=False)

print("Finished!")
print(f"Per-read file saved to: {PER_READ_OUT}")
print(f"Summary file saved to: {SUMMARY_OUT}")

