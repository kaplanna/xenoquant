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

context_window = 20  # ±10 bases

# --- Load XNA positions independently ---
mod_xna_pos = int(pd.read_csv(mod_bed, sep="\t", header=None).iloc[0, 1])
can_xna_pos = int(pd.read_csv(can_bed, sep="\t", header=None).iloc[0, 1])
print(f"📍 MOD XNA position: {mod_xna_pos}")
print(f"📍 CAN XNA position: {can_xna_pos}")

# --- Helper: extract signal window centered on XNA ---
def extract_signal_chunks(df, xna_pos, flank):
    signal_chunks = []

    for _, row in df.iterrows():
        try:
            ref_to_signal = row["Ref to Signal"]
            norm_signal = row["Normalized Signal"]

            if isinstance(ref_to_signal, str):
                ref_to_signal = ast.literal_eval(ref_to_signal)
            if isinstance(norm_signal, str):
                norm_signal = ast.literal_eval(norm_signal)

            ref_to_signal = np.array(ref_to_signal)
            norm_signal = np.array(norm_signal)


            # Ensure window fits in reference
            if xna_pos - flank < 0 or xna_pos + flank >= len(ref_to_signal):
                continue

            window = ref_to_signal[xna_pos - flank : xna_pos + flank + 1]
            if any(i is None for i in window):
                continue

            signal_indices = list(map(int, window))
            signal_chunk = norm_signal[signal_indices]

            signal_chunks.append(signal_chunk)

        except Exception as e:
            print(f"⚠️ Skipping read due to error: {e}")
            continue

    return np.vstack(signal_chunks) if signal_chunks else np.empty((0, 2 * flank + 1))


# --- Load CSVs ---
print("📄 Loading mod and can CSVs...")
mod_df = pd.read_csv(mod_csv)
can_df = pd.read_csv(can_csv)

# --- Extract aligned signal chunks ---
print("📦 Extracting signal chunks for mod...")
mod_chunks = extract_signal_chunks(mod_df, mod_xna_pos, context_window)

print("📦 Extracting signal chunks for can...")
can_chunks = extract_signal_chunks(can_df, can_xna_pos, context_window)

positions = np.arange(-context_window, context_window + 1)
plt.figure(figsize=(8, 4))

def plot_all_reads(data, color, label):
    for row in data:
        plt.plot(positions, row, color=color, alpha=0.1, linewidth=0.8)
    plt.plot([], [], color=color, alpha=0.6, label=label)  # dummy line for legend

plot_all_reads(can_chunks, "#1B9E77", "Canonical")     # Teal
plot_all_reads(mod_chunks, "#D95F02", "Modified")      # Magenta



plt.axvline(0, color='gray', linestyle='--', linewidth=1)
plt.xlabel("Position relative to XNA site")
plt.ylabel("Normalized signal")
plt.title(f"Per-read signal around {XNA} site (±{context_window} bp)")
plt.legend()
plt.tight_layout()
plt.show()

