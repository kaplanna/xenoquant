import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ast

# --- Set working directory and base type ---
XNA = 'S'

mod_working_dir = '/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/Manuscript_Model_Testing/DsPx/GDsA/Mod/PxT-Basecall'
can_working_dir = '/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/Manuscript_Model_Testing/DsPx/GDsA/Can/PxT-Basecall'

mod_csv = os.path.join(mod_working_dir, 'vis_extract', 'vis_output.csv')
can_csv = os.path.join(can_working_dir, 'vis_extract', 'vis_output.csv')

mod_bed = os.path.join(mod_working_dir, 'references', f'{XNA}.bed')
can_bed = os.path.join(can_working_dir, 'references', f'{XNA}.bed')

flank = 10  # ±10 bases

# --- Load XNA positions from BED files ---
mod_xna = int(pd.read_csv(mod_bed, sep="\t", header=None).iloc[0, 1])
can_xna = int(pd.read_csv(can_bed, sep="\t", header=None).iloc[0, 1])

print(f"📍 MOD XNA pos: {mod_xna} | CAN XNA pos: {can_xna}")

# --- Signal extraction ---
def extract_dacs_around_xna(df, xna_pos, flank):
    chunks = []
    for _, row in df.iterrows():
        try:
            dacs = np.array(list(map(int, row["Signal"].split(","))))
            ref_to_signal = row["Ref to Signal"]
            if isinstance(ref_to_signal, str):
                ref_to_signal = ast.literal_eval(ref_to_signal)
            ref_to_signal = np.array(ref_to_signal)

            if xna_pos - flank < 0 or xna_pos + flank >= len(ref_to_signal):
                continue

            sig_idx_window = ref_to_signal[xna_pos - flank : xna_pos + flank + 1]
            if any(x is None for x in sig_idx_window):
                continue

            signal_chunk = dacs[list(map(int, sig_idx_window))]
            chunks.append(signal_chunk)

        except Exception as e:
            print(f"⚠️ Skipped a row: {e}")
            continue

    return np.vstack(chunks) if chunks else np.empty((0, 2 * flank + 1))

# --- Load and extract ---
mod_df = pd.read_csv(mod_csv)
can_df = pd.read_csv(can_csv)

mod_chunks = extract_dacs_around_xna(mod_df, mod_xna, flank)
can_chunks = extract_dacs_around_xna(can_df, can_xna, flank)

positions = np.arange(-flank, flank + 1)

# =======================
# 🔹 PLOT 1: SPAGHETTI
# =======================
plt.figure(figsize=(9, 4))

def plot_spaghetti(chunks, color, label):
    for row in chunks:
        plt.plot(positions, row, color=color, alpha=0.1, linewidth=0.8)
    plt.plot([], [], color=color, alpha=0.6, label=label)

plot_spaghetti(can_chunks, "#4E79A7", "Canonical")
plot_spaghetti(mod_chunks, "#E15759", "Modified")

plt.axvline(0, color='gray', linestyle='--')
plt.xlabel("Position relative to XNA site")
plt.ylabel("Raw DAC signal")
plt.title("Raw signal (DAC) around XNA site — All Reads")
plt.legend()
plt.tight_layout()
plt.show()

# =======================
# 🔹 PLOT 2: MEAN ± STD
# =======================
plt.figure(figsize=(9, 4))

def plot_mean_std(chunks, color, label):
    if chunks.size == 0:
        return
    mean = np.mean(chunks, axis=0)
    std = np.std(chunks, axis=0)
    plt.plot(positions, mean, label=label, color=color)
    plt.fill_between(positions, mean - std, mean + std, color=color, alpha=0.3)

plot_mean_std(can_chunks, "#4E79A7", "Canonical")
plot_mean_std(mod_chunks, "#E15759", "Modified")

plt.axvline(0, color='gray', linestyle='--')
plt.xlabel("Position relative to XNA site")
plt.ylabel("Raw DAC signal")
plt.title("Raw signal (DAC) around XNA site — Mean ± Std")
plt.legend()
plt.tight_layout()
plt.show()

