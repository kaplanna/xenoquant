#!/usr/bin/env python3

import pandas as pd
import os

# ============================
# USER SETTINGS
# ============================
INPUT_TSV = "/home/marchandlab/DataAnalysis/Kaplan/basecall/8letter/251013_P8_H4_Basecall/PG-75_75_Basecall/remora_outputs/per-read_modifications.tsv"
OUTPUT_TSV = "/home/marchandlab/DataAnalysis/Kaplan/basecall/8letter/251013_P8_H4_Basecall/PG-75_75_Basecall/remora_outputs/per-read_modifications_focus-Filtered.tsv"

FOCUS_MIN = 40   # inclusive
FOCUS_MAX = 60   # inclusive

# ============================
# Load TSV
# ============================
df = pd.read_csv(INPUT_TSV, sep="\t")

# ============================
# Sanity check
# ============================
if "read_focus_base" not in df.columns:
    raise ValueError("Column 'read_focus_base' not found in TSV")

# ============================
# Filter by focus position
# ============================
filtered_df = df[
    (df["read_focus_base"] >= FOCUS_MIN) &
    (df["read_focus_base"] <= FOCUS_MAX)
].copy()

# ============================
# Save filtered TSV
# ============================
filtered_df.to_csv(OUTPUT_TSV, sep="\t", index=False)

print("Focus-position filtering complete")
print(f"Input reads:  {len(df)}")
print(f"Output reads: {len(filtered_df)}")
print(f"Saved to:     {OUTPUT_TSV}")

