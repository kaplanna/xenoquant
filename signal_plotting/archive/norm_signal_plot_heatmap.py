import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ast

# --- Config ---
XNA = 'B'
context_window = 1000
flank = context_window // 2
DOWNSAMPLE_FACTOR = 2

# --- Working directories ---
mod_working_dir = '/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/Manuscript_Model_Testing/BS/CBG/Mod/BA-Basecall'
can_working_dir = '/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/Manuscript_Model_Testing/BS/CBG/Can/BA-Basecall'

mod_csv = os.path.join(mod_working_dir, 'vis_extract', 'vis_output.csv')
can_csv = os.path.join(can_working_dir, 'vis_extract', 'vis_output.csv')

mod_bed = os.path.join(mod_working_dir, 'references', f'{XNA}.bed')
can_bed = os.path.join(can_working_dir, 'references', f'{XNA}.bed')

# --- Load XNA positions from BED files ---
mod_xna_pos = int(pd.read_csv(mod_bed, sep="\t", header=None).iloc[0, 1])
can_xna_pos = int(pd.read_csv(can_bed, sep="\t", header=None).iloc[0, 1])
print(f"📍 MOD XNA position: {mod_xna_pos}")
print(f"📍 CAN XNA position: {can_xna_pos}")

# --- Output path ---
vis_dir = os.path.join(mod_working_dir, 'vis_extract')
os.makedirs(vis_dir, exist_ok=True)
fig_path = os.path.join(vis_dir, f"{XNA}_norm_signal_heatmap_side_by_side.pdf")

# --- Downsample helper ---
def downsample(matrix, factor=10):
    return np.array([
        row[:len(row) // factor * factor].reshape(-1, factor).mean(axis=1)
        for row in matrix
    ]) if matrix.size else np.empty((0,))

# --- Extract normalized signal chunks around BED-defined XNA ---
def extract_norm_chunks(df, xna_pos, flank):
    chunks = []
    for _, row in df.iterrows():
        try:
            norm_signal = row["Normalized Signal"]
            if isinstance(norm_signal, str):
                norm_signal = ast.literal_eval(norm_signal)
            norm_signal = np.array(norm_signal)

            ref_to_signal = row["Ref to Signal"]
            if isinstance(ref_to_signal, str):
                ref_to_signal = ast.literal_eval(ref_to_signal)
            ref_to_signal = np.array(ref_to_signal)

            sig_center = ref_to_signal[xna_pos]
            if sig_center is None or sig_center - flank < 0 or sig_center + flank >= len(norm_signal):
                continue

            chunk = norm_signal[sig_center - flank : sig_center + flank + 1]
            chunks.append(chunk)

        except Exception as e:
            print(f"⚠️ Skipping read due to: {e}")
            continue

    return np.vstack(chunks) if chunks else np.empty((0, 2 * flank + 1))

# --- Load data and extract signal chunks ---
mod_df = pd.read_csv(mod_csv)
can_df = pd.read_csv(can_csv)

mod_chunks = extract_norm_chunks(mod_df, mod_xna_pos, flank)
can_chunks = extract_norm_chunks(can_df, can_xna_pos, flank)

# --- Sanity check ---
if mod_chunks.size == 0 or can_chunks.size == 0:
    print("⚠️ No signal chunks extracted. Exiting.")
    exit()

# --- Downsample ---
mod_chunks_ds = downsample(mod_chunks, DOWNSAMPLE_FACTOR)
can_chunks_ds = downsample(can_chunks, DOWNSAMPLE_FACTOR)
positions_ds = np.linspace(-flank, flank, mod_chunks_ds.shape[1])

# --- Side-by-side heatmap plot ---
fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

# --- Modified ---
im0 = axes[0].imshow(mod_chunks_ds, aspect='auto', cmap='magma',
                     extent=[-flank, flank, 0, mod_chunks_ds.shape[0]])
axes[0].axvline(0, color='white', linestyle='--')
axes[0].set_title("Modified")
axes[0].set_xlabel("Signal index (downsampled)")
axes[0].set_ylabel("Read index")

# --- Canonical ---
im1 = axes[1].imshow(can_chunks_ds, aspect='auto', cmap='magma',
                     extent=[-flank, flank, 0, can_chunks_ds.shape[0]])
axes[1].axvline(0, color='white', linestyle='--')
axes[1].set_title("Canonical")
axes[1].set_xlabel("Signal index (downsampled)")

# --- Colorbar ---
cbar = fig.colorbar(im1, ax=axes.ravel().tolist(), shrink=0.85, label="Normalized signal")

# --- Save and show ---
plt.suptitle(f"Normalized signal around XNA site ({XNA})", fontsize=14)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig(fig_path, bbox_inches='tight')
plt.show()

print(f"✅ Side-by-side heatmap saved to: {fig_path}")

