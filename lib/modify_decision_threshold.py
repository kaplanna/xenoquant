#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os

# ============================================================
# USER SETTINGS
# ============================================================

DEMUX_DIR = "/home/marchandlab/DataAnalysis/Kaplan/basecall/8letter/251214_H4_BC/ZC-100-100_full/demux"   # <-- EDIT
INPUT_TSV = os.path.join(DEMUX_DIR, "demux_per-read_modifications.tsv")

THRESHOLD = 0.95   # must be > 0.5

# Outputs
PER_READ_OUT = os.path.join(
    DEMUX_DIR,
    f"demux_per-read_modifications_recalled_T{THRESHOLD:.2f}.tsv"
)

OVERALL_OUT = os.path.join(
    DEMUX_DIR,
    f"overall_demux_results_recalled_T{THRESHOLD:.2f}.csv"
)

# ============================================================
# Helpers
# ============================================================

def parse_p1(prob_string):
    """Extract P(class1) from 'p0,p1'"""
    try:
        return float(prob_string.split(",")[1])
    except Exception:
        return np.nan


# ============================================================
# Main
# ============================================================

print(f"[INFO] Loading {INPUT_TSV}")
df = pd.read_csv(INPUT_TSV, sep="\t")

# Extract P(class1)
df["p_class1"] = df["class_probs"].apply(parse_p1)

# HARD reassignment (no ambiguous state)
df["class_pred_reassigned"] = (df["p_class1"] >= THRESHOLD).astype(int)

# ------------------------------------------------------------
# Save per-read reassigned output
# ------------------------------------------------------------
df.to_csv(PER_READ_OUT, sep="\t", index=False)

print("[INFO] Per-read reassigned file written:")
print(f"       {PER_READ_OUT}")

# ------------------------------------------------------------
# Recompute overall results
# ------------------------------------------------------------
overall = (
    df
    .groupby(["sample_id", "barcode_pair"])
    .agg(
        Total_Reads=("read_id", "size"),
        Number_of_1s=("class_pred_reassigned", lambda x: (x == 1).sum()),
        Number_of_0s=("class_pred_reassigned", lambda x: (x == 0).sum())
    )
    .reset_index()
)

overall["Percentage_1"] = 100 * overall["Number_of_1s"] / overall["Total_Reads"]
overall["Percentage_0"] = 100 * overall["Number_of_0s"] / overall["Total_Reads"]

overall.to_csv(OVERALL_OUT, index=False)

print("[DONE]")
print(f"  Decision threshold: {THRESHOLD}")
print(f"  Overall results written to:")
print(f"    {OVERALL_OUT}")

# ------------------------------------------------------------
# Summary stats
# ------------------------------------------------------------
changed = (df["class_pred"] != df["class_pred_reassigned"]).sum()
print("\n[SUMMARY]")
print(f"  Reads changed vs Remora argmax: {changed} / {len(df)}")

