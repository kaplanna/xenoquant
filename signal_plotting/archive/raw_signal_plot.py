import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ast

# --- Set working directory and base type ---
XNA = 'Z'

mod_working_dir = '/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/Manuscript_Model_Testing/PZ/Mod/ZC-Basecall'
can_working_dir = '/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/Manuscript_Model_Testing/PZ/Can/ZC-Basecall'

mod_csv = os.path.join(mod_working_dir, 'vis_extract', 'vis_output.csv')
can_csv = os.path.join(can_working_dir, 'vis_extract', 'vis_output.csv')


mod_bed = os.path.join(mod_working_dir, 'references', f'{XNA}.bed')
can_bed = os.path.join(can_working_dir, 'references', f'{XNA}.bed')

context_window = 500  # ± signal samples
flank = context_window // 2

# Output figure path
vis_dir = os.path.join(mod_working_dir, 'vis_extract')  # or can_working_dir
os.makedirs(vis_dir, exist_ok=True)

spaghetti_path = os.path.join(vis_dir, f"{XNA}_raw_signal_spaghetti.pdf")
meanstd_path = os.path.join(vis_dir, f"{XNA}_raw_signal_mean_std.pdf")


# --- Load XNA positions from BED files ---
mod_xna_pos = int(pd.read_csv(mod_bed, sep="\t", header=None).iloc[0, 1])
can_xna_pos = int(pd.read_csv(can_bed, sep="\t", header=None).iloc[0, 1])
print(f"📍 MOD XNA position: {mod_xna_pos}")
print(f"📍 CAN XNA position: {can_xna_pos}")

def downsample(matrix, factor=10):
    # Reduce each row to 1/factor length by averaging chunks
    return np.array([
        row[:len(row) // factor * factor].reshape(-1, factor).mean(axis=1)
        for row in matrix
    ])


# --- Extract normalized signal window around BED-defined XNA ---
def extract_dac_signal_chunks(df, xna_pos, flank):
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
            print(f"⚠️ Skipping row due to: {e}")
            continue

    return np.vstack(chunks) if chunks else np.empty((0, 2 * flank + 1))


def plot_spaghetti(chunks, color, label):
    for row in chunks:
        plt.plot(positions_ds, row, color=color, alpha=0.3, linewidth=0.8)
    plt.plot([], [], color=color, alpha=0.6, label=label)

def plot_mean_std(chunks, color, label):
    if chunks.size == 0:
        return
    mean = np.mean(chunks, axis=0)
    std = np.std(chunks, axis=0)
    plt.plot(positions_ds, mean, label=label, color=color)
    plt.fill_between(positions_ds, mean - std, mean + std, color=color, alpha=0.3)

# --- Load datasets ---
mod_df = pd.read_csv(mod_csv)
can_df = pd.read_csv(can_csv)

mod_chunks = extract_dac_signal_chunks(mod_df, mod_xna_pos, flank)
can_chunks = extract_dac_signal_chunks(can_df, can_xna_pos, flank)


# Downsample
DOWNSAMPLE_FACTOR = 5
mod_chunks_ds = downsample(mod_chunks, factor=DOWNSAMPLE_FACTOR)
can_chunks_ds = downsample(can_chunks, factor=DOWNSAMPLE_FACTOR)

positions_ds = np.linspace(-flank, flank, mod_chunks_ds.shape[1])

# Spaghetti plot
plt.figure(figsize=(10, 4))
plot_spaghetti(can_chunks_ds, "#4E79A7", "Canonical")
plot_spaghetti(mod_chunks_ds, "#E15759", "Modified")
plt.axvline(0, color='gray', linestyle='--')
plt.xlabel("Signal index (downsampled)")
plt.ylabel("Raw DAC signal")
plt.title(f"Raw DAC signal — All Reads ({XNA})")
plt.legend()
plt.tight_layout()
plt.savefig(spaghetti_path, bbox_inches='tight')
plt.show()


# Mean ± std plot
plt.figure(figsize=(10, 4))
plot_mean_std(can_chunks_ds, "#4E79A7", "Canonical")
plot_mean_std(mod_chunks_ds, "#E15759", "Modified")
plt.axvline(0, color='gray', linestyle='--')
plt.xlabel("Signal index (downsampled)")
plt.ylabel("Raw DAC signal")
plt.title(f"Raw DAC signal — Mean ± Std ({XNA})")
plt.legend()
plt.tight_layout()
plt.savefig(meanstd_path, bbox_inches='tight')
plt.show()


