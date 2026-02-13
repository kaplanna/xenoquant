#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os

# ============================
# User settings
# ============================
PER_READ_CSV = "/home/marchandlab/DataAnalysis/Kaplan/basecall/8letter/251013_P8_H4_Basecall/251213_H4-SZ-Combined_Calc/per_read_comparison.csv"
OUTDIR = "/home/marchandlab/DataAnalysis/Kaplan/basecall/8letter/251013_P8_H4_Basecall/251213_H4-SZ-Combined_Calc/plots"
PLOT_NAME = "focus_position_histogram"

BARCODE_FILTER_CSV = '/home/marchandlab/DataAnalysis/Kaplan/basecall/8letter/251013_P8_H4_Basecall/251213_H4-PB-Combined_Calc/H4_barcodes.csv'

BINS = 80

os.makedirs(OUTDIR, exist_ok=True)

# ============================
# Matplotlib styling
# ============================
mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 11,
    "pdf.fonttype": 42,
    "ps.fonttype": 42
})

# ============================
# Load data
# ============================
df = pd.read_csv(PER_READ_CSV)

# ============================
# Optional barcode filtering
# ============================
if BARCODE_FILTER_CSV is not None:
    bc_df = pd.read_csv(BARCODE_FILTER_CSV, header=None)
    barcode_col = bc_df.columns[0]
    barcodes_to_keep = set(bc_df[barcode_col].astype(str))

    n_before = len(df)
    df = df[df["barcode_pair"].astype(str).isin(barcodes_to_keep)]
    n_after = len(df)

    print(f"Barcode filtering enabled:")
    print(f"  Reads before: {n_before}")
    print(f"  Reads after:  {n_after}")

# ============================
# Detect focus column names
# ============================
if "focus_base_m1" in df.columns:
    f1, f2 = "focus_base_m1", "focus_base_m2"
elif "read_focus_base_m1" in df.columns:
    f1, f2 = "read_focus_base_m1", "read_focus_base_m2"
else:
    raise ValueError("Focus base columns not found.")

df = df.dropna(subset=[f1, f2])

focus_m1 = df[f1].astype(int)
focus_m2 = df[f2].astype(int)

# ============================
# Shared bins and limits
# ============================
min_pos = min(focus_m1.min(), focus_m2.min())
max_pos = max(focus_m1.max(), focus_m2.max())
bins = np.linspace(min_pos, max_pos, BINS)

max_count = max(
    np.histogram(focus_m1, bins=bins)[0].max(),
    np.histogram(focus_m2, bins=bins)[0].max()
)

# ============================
# Plot
# ============================
fig, axes = plt.subplots(
    nrows=2, ncols=1,
    figsize=(4.5, 5.5),
    sharex=True,
    sharey=True
)

# ---- Model 1 ----
axes[0].hist(
    focus_m1,
    bins=bins,
    edgecolor="black",
    linewidth=0.6
)
axes[0].set_title("Model 1", pad=6)
axes[0].set_ylabel("Read count")
axes[0].set_ylim(0, max_count * 1.05)

# ---- Model 2 ----
axes[1].hist(
    focus_m2,
    bins=bins,
    edgecolor="black",
    linewidth=0.6
)
axes[1].set_title("Model 2", pad=6)
axes[1].set_ylabel("Read count")
axes[1].set_xlabel("Focus base position (basecalled sequence)")
axes[1].set_ylim(0, max_count * 1.05)

fig.tight_layout()

# ============================
# Save (PDF only)
# ============================
out_pdf = os.path.join(OUTDIR, f"{PLOT_NAME}.pdf")
fig.savefig(out_pdf, bbox_inches="tight")
plt.close(fig)

print("Saved:")
print(out_pdf)

