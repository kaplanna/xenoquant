import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ast

# --- Set working directory and base type ---
XNA = 'S'

mod_working_dir = '/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/Manuscript_Model_Testing/BS/CBG/Mod/ST-Basecall'
can_working_dir = '/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/Manuscript_Model_Testing/BS/CBG/Can/ST-Basecall'

mod_csv = os.path.join(mod_working_dir, 'vis_extract', 'vis_output.csv')
can_csv = os.path.join(can_working_dir, 'vis_extract', 'vis_output.csv')

mod_bed = os.path.join(mod_working_dir, 'references', f'{XNA}.bed')
can_bed = os.path.join(can_working_dir, 'references', f'{XNA}.bed')

context_window = 500
flank = context_window // 2

# --- Load XNA positions from BED files ---
mod_xna_pos = int(pd.read_csv(mod_bed, sep="\t", header=None).iloc[0, 1])
can_xna_pos = int(pd.read_csv(can_bed, sep="\t", header=None).iloc[0, 1])
print(f"📍 MOD XNA position: {mod_xna_pos}")
print(f"📍 CAN XNA position: {can_xna_pos}")

# --- Downsample helper ---
def downsample(matrix, factor=5):
    return np.array([
        row[:len(row) // factor * factor].reshape(-1, factor).mean(axis=1)
        for row in matrix
    ]) if matrix.size else np.empty((0,))

# --- Extract normalized signal chunks ---
def extract_norm_signal_chunks(df, xna_pos, flank):
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
            print(f"⚠️ Skipping row due to: {e}")
            continue

    return np.vstack(chunks) if chunks else np.empty((0, 2 * flank + 1))

# --- Plotting helpers ---
def plot_spaghetti(chunks, color, label):
    for row in chunks:
        plt.plot(positions_ds, row, color=color, alpha=0.15, linewidth=0.8)
    plt.plot([], [], color=color, alpha=0.8, label=label)


def plot_median_iqr(chunks, color, fill_color, label):
    if chunks.size == 0:
        return
    median = np.median(chunks, axis=0)
    q25 = np.percentile(chunks, 25, axis=0)
    q75 = np.percentile(chunks, 75, axis=0)
    plt.plot(positions_ds, median, label=label, color=color, linewidth=2)
    plt.fill_between(positions_ds, q25, q75, color=fill_color, alpha=0.3)

# --- Load and process data ---
mod_df = pd.read_csv(mod_csv)
can_df = pd.read_csv(can_csv)

mod_chunks = extract_norm_signal_chunks(mod_df, mod_xna_pos, flank)
can_chunks = extract_norm_signal_chunks(can_df, can_xna_pos, flank)

DOWNSAMPLE_FACTOR = 2
mod_chunks_ds = downsample(mod_chunks, factor=DOWNSAMPLE_FACTOR)
can_chunks_ds = downsample(can_chunks, factor=DOWNSAMPLE_FACTOR)

positions_ds = np.linspace(-flank, flank, mod_chunks_ds.shape[1])

# --- Output figure paths ---
vis_dir = os.path.join(mod_working_dir, 'vis_extract')
os.makedirs(vis_dir, exist_ok=True)

spaghetti_path = os.path.join(vis_dir, f"{XNA}_norm_signal_spaghetti.pdf")
iqr_path = os.path.join(vis_dir, f"{XNA}_norm_signal_median_iqr.pdf")

# --- Plot 1: Spaghetti ---
plt.figure(figsize=(10, 4))
plot_spaghetti(can_chunks_ds, "#000000", "Canonical")
plot_spaghetti(mod_chunks_ds, "#D00072", "Modified")
plt.axvline(0, color='gray', linestyle='--')
plt.xlabel("Signal index (downsampled)")
plt.ylabel("Normalized signal")
plt.title("Normalized signal — All Reads")
plt.legend()
plt.tight_layout()
plt.savefig(spaghetti_path, bbox_inches='tight')
plt.show()


plt.figure(figsize=(10, 4))
plot_median_iqr(can_chunks_ds, "#353535", "#999999", "Canonical")
plot_median_iqr(mod_chunks_ds, "#D00072", "#CE6AA6", "Modified")
plt.axvline(0, color='gray', linestyle='--')
plt.xlabel("Signal index (downsampled)")
plt.ylabel("Normalized signal")
plt.title("Normalized signal — Median ± IQR")
plt.legend()
plt.tight_layout()
plt.savefig(iqr_path, bbox_inches='tight')
plt.show()

