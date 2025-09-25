import pandas as pd
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, roc_auc_score
from xr_params import *



def define_directories(base_path):
    results_dir = os.path.join(base_path, 'alignment_results')
    bed_dir = os.path.join(base_path, 'references')
    os.makedirs(results_dir, exist_ok=True)
    return results_dir, bed_dir

def read_bed_file(bed_dir):
    bed_path = os.path.join(bed_dir, f'{mod_base}.bed')
    if not os.path.exists(bed_path):
        print(f"Missing BED file at: {bed_path}")
        sys.exit(1)
    df = pd.read_csv(bed_path, delim_whitespace=True, header=None)
    df.columns = ['Alignment', 'Start', 'End', 'Name', 'Score', 'Strand']
    return df

def read_modifications_file(base_path):
    mod_file = os.path.join(base_path, 'remora_outputs/per-read_modifications.tsv')
    if not os.path.exists(mod_file):
        print(f"Missing modifications file at: {mod_file}")
        sys.exit(1)
    df = pd.read_csv(mod_file, sep='\t')
    df.columns = ['read_id', 'read_focus_base', 'label', 'class_pred', 'class_probs',
                  'reference_sequence', 'flag', 'ref_start_pos', 'cigar_string',
                  'ref_length', 'basecalled_sequence', 'q_score']
    return df

def process_data(df, alignment):
    df = df[df['reference_sequence'] == alignment].copy()
    if df.empty:
        return df
    df[['class_0_probs', 'class_1_probs']] = df['class_probs'].str.split(',', expand=True).astype(float)
    df['class_pred'] = pd.to_numeric(df['class_pred'], errors='coerce').fillna(0).astype(int)
    return df[['read_id', 'class_pred', 'class_0_probs', 'class_1_probs']].drop_duplicates('read_id')

def calculate_results(df, alignment):
    n1 = df['class_pred'].sum()
    n0 = len(df) - n1
    return {
        'Sequence': alignment,
        'Total Alignments': len(df),
        'Number of 1s': n1,
        'Number of 0s': n0,
        'Percentage 1': round(100 * n1 / len(df), 2),
        'Percentage 0': round(100 * n0 / len(df), 2)
    }

def plot_class_probs(df, title, save_path, bins=100):
    if df.empty:
        return
    plt.figure(figsize=(8, 5))
    plt.hist(df['class_1_probs'], bins=bins, color='steelblue', edgecolor='black')
    plt.xlabel('Class 1 Probability')
    plt.ylabel('Read Count')
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()

def save_alignment_df(df, alignment, results_dir):
    filename = f"alignment_results_{alignment}.csv"
    path = os.path.join(results_dir, filename)
    df.to_csv(path, index=False)
    return filename

def analyze_results(base_path):
    results_dir, bed_dir = define_directories(base_path)
    df_bed = read_bed_file(bed_dir)
    df_mods = read_modifications_file(base_path)

    full_results = []
    full_alignment_results_path = os.path.join(results_dir, "full_alignment_results.csv")
    if os.path.exists(full_alignment_results_path):
        os.remove(full_alignment_results_path)

    for _, row in df_bed.iterrows():
        alignment = row["Alignment"]
        df = process_data(df_mods, alignment)
        if df.empty:
            print(f"Skipping {alignment} (no reads)")
            continue

        # Save individual alignment
        df.to_csv(full_alignment_results_path, mode='a',
                  header=not os.path.exists(full_alignment_results_path),
                  index=False)

        save_alignment_df(df, alignment, results_dir)
        result = calculate_results(df, alignment)
        full_results.append(result)

        # Plot
        plot_path = os.path.join(results_dir, f'class1_probability_distribution_{alignment}.png')
        plot_class_probs(df, f"Class 1 Probability Distribution: {alignment}", plot_path)
        print(f"Saved: {plot_path}")

    # Save summary
    summary_df = pd.DataFrame(full_results)
    summary_df.to_csv(os.path.join(results_dir, 'overall_results.csv'), index=False)
    print("Analysis complete. Results saved to alignment_results/")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python xr_results.py /path/to/base_directory")
        sys.exit(1)

    base_path = os.path.expanduser(sys.argv[1])
    analyze_results(base_path)

