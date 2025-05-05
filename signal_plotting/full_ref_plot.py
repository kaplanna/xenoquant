import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ast

# --- User Configuration ---
XNA = "B"                         # Your modified base
FLANK = 10                        # Bases to highlight around the XNA
MAX_READS = 100                    # Max reads to plot per dataset

mod_working_dir = '/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/Manuscript_Model_Testing/BS/CBG/Mod/BA-Basecall'
can_working_dir = '/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/Manuscript_Model_Testing/BS/CBG/Can/BA-Basecall'

# --- Output figure paths ---
vis_dir = os.path.join(mod_working_dir, 'vis_extract')
os.makedirs(vis_dir, exist_ok=True)

mod_csv = os.path.join(mod_working_dir, 'vis_extract', 'vis_output.csv')
can_csv = os.path.join(can_working_dir, 'vis_extract', 'vis_output.csv')

mod_bed = os.path.join(mod_working_dir, 'references', f'{XNA}.bed')
can_bed = os.path.join(can_working_dir, 'references', f'{XNA}.bed')

# --- Load BED files ---
mod_xna_pos = int(pd.read_csv(mod_bed, sep="\t", header=None).iloc[0, 1])
can_xna_pos = int(pd.read_csv(can_bed, sep="\t", header=None).iloc[0, 1])

# --- Load CSVs ---
mod_df = pd.read_csv(mod_csv)
can_df = pd.read_csv(can_csv)

# --- Extract zoomed-in signal matrix ---
def extract_signal_matrix(df, xna_pos, flank):
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

            start = xna_pos - flank
            end = xna_pos + flank + 1

            if start < 0 or end > len(ref_to_signal):
                continue

            window = ref_to_signal[start:end]
            if any(i is None for i in window):
                continue

            signal_indices = [int(i) for i in window]
            if max(signal_indices) >= len(norm_signal):
                continue

            chunk = norm_signal[signal_indices]
            signal_chunks.append(chunk)

        except Exception as e:
            print(f"⚠️ Skipping read: {e}")
            continue

    return np.vstack(signal_chunks) if signal_chunks else np.empty((0, 2 * flank + 1))

# --- Full-length spaghetti plot ---
def plot_full_ref_signals(df, xna_pos, label, save_path):
    plt.figure(figsize=(12, 4))
    plotted = 0

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

            signal = []
            ref_pos = []

            for i, idx in enumerate(ref_to_signal):
                if idx is None or idx >= len(norm_signal):
                    continue
                signal.append(norm_signal[idx])
                ref_pos.append(i)

            if not signal:
                continue

            plt.plot(ref_pos, signal, alpha=0.2, linewidth=0.6, color='gray')
            plotted += 1
            if plotted >= MAX_READS:
                break

        except Exception as e:
            print(f"⚠️ Skipping full-ref read: {e}")
            continue

    plt.axvline(xna_pos, linestyle='--', color='black', label='XNA site')
    plt.axvspan(xna_pos - FLANK, xna_pos + FLANK, color='gray', alpha=0.1, label='±10 bp')
    plt.title(f"{label} — Full Reference-Aligned Signal")
    plt.xlabel("Reference Position (local)")
    plt.ylabel("Normalized Signal")
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path)
    plt.show()

# --- Zoomed spaghetti + mean ± std ---
def plot_overlay_zoomed_signals(mod_data, can_data, label, save_path):
    positions = np.arange(-FLANK, FLANK + 1)
    plt.figure(figsize=(8, 4))

    # Means + stds
    mod_mean = np.mean(mod_data, axis=0)
    mod_std = np.std(mod_data, axis=0)
    can_mean = np.mean(can_data, axis=0)
    can_std = np.std(can_data, axis=0)

    # Canonical (green)
    plt.plot(positions, can_mean, color="#1B9E77", label="Canonical")
    plt.fill_between(positions, can_mean - can_std, can_mean + can_std, color="#1B9E77", alpha=0.2)

    # Modified (orange)
    plt.plot(positions, mod_mean, color="#D95F02", label="Modified")
    plt.fill_between(positions, mod_mean - mod_std, mod_mean + mod_std, color="#D95F02", alpha=0.2)

    plt.axvline(0, linestyle='--', color='gray')
    plt.title(f"{label} — Zoomed Signal ±{FLANK} bp")
    plt.xlabel("Position relative to XNA site")
    plt.ylabel("Normalized Signal")
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path)
    plt.show()
def plot_overlay_step_signal(mod_data, can_data, label, save_path):
    positions = np.arange(-FLANK, FLANK + 1)
    mod_mean = np.mean(mod_data, axis=0)
    can_mean = np.mean(can_data, axis=0)

    step_x = []
    mod_y = []
    can_y = []

    for i in range(len(positions) - 1):
        x0, x1 = positions[i], positions[i + 1]
        step_x.extend([x0, x1])
        mod_y.extend([mod_mean[i], mod_mean[i]])
        can_y.extend([can_mean[i], can_mean[i]])

    plt.figure(figsize=(8, 4))

    # Plot step lines
    plt.plot(step_x, mod_y, color="#D95F02", label="Modified (Step)")
    plt.plot(step_x, can_y, color="#1B9E77", label="Canonical (Step)")

    plt.axvline(0, linestyle='--', color='gray')
    plt.title(f"{label} — Stepwise Average Signal ±{FLANK} bp")
    plt.xlabel("Position relative to XNA site")
    plt.ylabel("Normalized Signal")
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path)
    plt.show()

def plot_signal_difference(mod_data, can_data, label, save_path):
    positions = np.arange(-FLANK, FLANK + 1)
    mod_mean = np.mean(mod_data, axis=0)
    can_mean = np.mean(can_data, axis=0)

    diff = mod_mean - can_mean

    plt.figure(figsize=(8, 4))
    plt.plot(positions, diff, color="#7570B3", label='Modified - Canonical')
    plt.axhline(0, color='gray', linestyle='--')
    plt.axvline(0, linestyle='--', color='gray')
    plt.title(f"{label} — Signal Difference (Mod - Can)")
    plt.xlabel("Position relative to XNA site")
    plt.ylabel("Δ Normalized Signal")
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path)
    plt.show()


# --- Process + Plot ---
print("📦 Extracting zoomed-in signals...")
mod_chunks = extract_signal_matrix(mod_df, mod_xna_pos, FLANK)
can_chunks = extract_signal_matrix(can_df, can_xna_pos, FLANK)

print("📈 Plotting full-length signals...")
plot_full_ref_signals(mod_df, mod_xna_pos, "Modified", os.path.join(vis_dir, "mod_full_signal.png"))
plot_full_ref_signals(can_df, can_xna_pos, "Canonical", os.path.join(vis_dir, "can_full_signal.png"))


print("📈 Plotting combined Canonical vs Modified overlay...")
plot_overlay_zoomed_signals(mod_chunks, can_chunks, "CBG ST", os.path.join(vis_dir, "overlay_zoom_signal.png"))

print("📈 Plotting stepwise Canonical vs Modified overlay...")
plot_overlay_step_signal(mod_chunks, can_chunks, "CBG ST", os.path.join(vis_dir, "overlay_step_signal.png"))

print("📈 Plotting signal difference (Modified - Canonical)...")
plot_signal_difference(mod_chunks, can_chunks, "CBG ST", os.path.join(vis_dir, "diff_signal.png"))


print("✅ All plots saved in:", vis_dir)

