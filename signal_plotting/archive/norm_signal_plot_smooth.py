import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ast
from scipy.ndimage import gaussian_filter1d



# --- Set working directory and base type ---
XNA = 'B'
context_window = 800  # total window size
flank = context_window // 2
DOWNSAMPLE_FACTOR = 2
SMOOTH_SIGMA = 5  # Gaussian smoothing level


mod_working_dir = '/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/Manuscript_Model_Testing/BS/CBG/Mod/BA-Basecall'
can_working_dir = '/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/Manuscript_Model_Testing/BS/CBG/Can/BA-Basecall'

mod_csv = os.path.join(mod_working_dir, 'vis_extract', 'vis_output.csv')
can_csv = os.path.join(can_working_dir, 'vis_extract', 'vis_output.csv')

mod_bed = os.path.join(mod_working_dir, 'references', f'{XNA}.bed')
can_bed = os.path.join(can_working_dir, 'references', f'{XNA}.bed')

context_window = 1000  # ±10 signal samples
flank = context_window // 2

# --- Load XNA positions from BED files ---
mod_xna_pos = int(pd.read_csv(mod_bed, sep="\t", header=None).iloc[0, 1])
can_xna_pos = int(pd.read_csv(can_bed, sep="\t", header=None).iloc[0, 1])
print(f"📍 MOD XNA position: {mod_xna_pos}")
print(f"📍 CAN XNA position: {can_xna_pos}")

# --- Output path ---
vis_dir = os.path.join(mod_working_dir, 'vis_extract')
os.makedirs(vis_dir, exist_ok=True)
fig_path = os.path.join(vis_dir, f"{XNA}_raw_signal_smoothed.pdf")

# --- Downsample helper ---
def downsample(matrix, factor=10):
    return np.array([
        row[:len(row) // factor * factor].reshape(-1, factor).mean(axis=1)
        for row in matrix
    ]) if matrix.size else np.empty((0,))

# --- Extract DAC signal chunks around BED-defined XNA ---
def extract_dac_chunks(df, xna_pos, flank):
    chunks = []
    for _, row in df.iterrows():
        try:
            dacs = row["Signal"]
            if isinstance(dacs, str):
                dacs = list(map(int, dacs.split(",")))
            dacs = np.array(dacs)

            ref_to_signal = row["Ref to Signal"]
            if isinstance(ref_to_signal, str):
                ref_to_signal = ast.literal_eval(ref_to_signal)
            ref_to_signal = np.array(ref_to_signal)

            sig_center = ref_to_signal[xna_pos]
            if sig_center is None or sig_center - flank < 0 or sig_center + flank >= len(dacs):
                continue

            chunk = dacs[sig_center - flank : sig_center + flank + 1]
            chunks.append(chunk)

        except Exception as e:
            print(f"⚠️ Skipping read due to: {e}")
            continue

    return np.vstack(chunks) if chunks else np.empty((0, 2 * flank + 1))

# --- Load data and extract signal chunks ---
mod_df = pd.read_csv(mod_csv)
can_df = pd.read_csv(can_csv)

mod_chunks = extract_dac_chunks(mod_df, mod_xna_pos, flank)
can_chunks = extract_dac_chunks(can_df, can_xna_pos, flank)

# --- Sanity check ---
if mod_chunks.size == 0 or can_chunks.size == 0:
    print("⚠️ No signal chunks extracted. Exiting.")
    exit()

# --- Downsample chunks ---
mod_chunks_ds = downsample(mod_chunks, DOWNSAMPLE_FACTOR)
can_chunks_ds = downsample(can_chunks, DOWNSAMPLE_FACTOR)
positions_ds = np.linspace(-flank, flank, mod_chunks_ds.shape[1])

# --- Smoothed mean plot ---
plt.figure(figsize=(10, 4))

mean_mod = np.mean(mod_chunks_ds, axis=0)
mean_can = np.mean(can_chunks_ds, axis=0)

smooth_mod = gaussian_filter1d(mean_mod, sigma=SMOOTH_SIGMA)
smooth_can = gaussian_filter1d(mean_can, sigma=SMOOTH_SIGMA)

plt.plot(positions_ds, smooth_can, label="Canonical", color="#4E79A7", linewidth=2)
plt.plot(positions_ds, smooth_mod, label="Modified", color="#E15759", linewidth=2)

plt.axvline(0, color='gray', linestyle='--')
plt.xlabel("Signal index (downsampled)")
plt.ylabel("Raw DAC signal")
plt.title(f"Smoothed mean raw signal around XNA site ({XNA})")
plt.legend()
plt.tight_layout()
plt.savefig(fig_path, bbox_inches='tight')
plt.show()

print(f"✅ Figure saved to: {fig_path}")

