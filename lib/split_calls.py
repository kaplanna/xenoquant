#!/usr/bin/env python3

import pandas as pd
import numpy as np

# =========================
# USER SETTINGS
# =========================

INPUT_FILE = "/home/marchandlab/DataAnalysis/Kaplan/basecall/8letter/251214_H4_BC/PP-125_125/demux/demux_per-read_modifications.tsv"   # TSV or CSV
OUTPUT_PER_READ = "/home/marchandlab/DataAnalysis/Kaplan/basecall/8letter/251214_H4_BC/PP-125_125/demux/per-read_with_calls.tsv"
OUTPUT_SUMMARY = "/home/marchandlab/DataAnalysis/Kaplan/basecall/8letter/251214_H4_BC/PP-125_125/demux/reference_sequence_call_summary.csv"
#!/usr/bin/env python3



SEP = "\t"

FOCUS_POS_1 = 57
FOCUS_POS_2 = 78

MAX_DISTANCE = None      # e.g. 6, or None
MIDPOINT_TOL = None      # e.g. 2, or None

# =========================
# Load data
# =========================

df = pd.read_csv(INPUT_FILE, sep=SEP)

required_cols = {
    "read_id",
    "read_focus_base",
    "class_pred",
    "barcode_pair"
}

missing = required_cols - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns: {missing}")

# =========================
# Distance calculations
# =========================

df["dist_1"] = (df["read_focus_base"] - FOCUS_POS_1).abs()
df["dist_2"] = (df["read_focus_base"] - FOCUS_POS_2).abs()

# Assign calls
df["call"] = np.where(df["dist_1"] <= df["dist_2"], "call_1", "call_2")

# =========================
# Optional filtering
# =========================

if MAX_DISTANCE is not None:
    df = df[(df["dist_1"] <= MAX_DISTANCE) | (df["dist_2"] <= MAX_DISTANCE)]

if MIDPOINT_TOL is not None:
    midpoint = (FOCUS_POS_1 + FOCUS_POS_2) / 2
    df = df[(df["read_focus_base"] - midpoint).abs() > MIDPOINT_TOL]

print(f"[INFO] Reads retained after filtering: {len(df)}")

# =========================
# Save per-read output
# =========================

df.to_csv(OUTPUT_PER_READ, sep="\t", index=False)
print(f"[OK] Saved per-read calls → {OUTPUT_PER_READ}")

# =========================
# Summary per barcode pair & call
# =========================

summary = (
    df.groupby(["barcode_pair", "call"])
      .agg(
          Total_Reads=("read_id", "size"),
          Number_of_1s=("class_pred", lambda x: (x == 1).sum()),
          Number_of_0s=("class_pred", lambda x: (x == 0).sum())
      )
      .reset_index()
)

summary["Percentage_1"] = 100 * summary["Number_of_1s"] / summary["Total_Reads"]
summary["Percentage_0"] = 100 * summary["Number_of_0s"] / summary["Total_Reads"]

summary.to_csv(OUTPUT_SUMMARY, index=False)
print(f"[OK] Saved summary → {OUTPUT_SUMMARY}")


