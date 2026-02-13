import pandas as pd
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, roc_auc_score
from xr_params import *

# ======================================
# USER SETTINGS — HARD-CODED DIRECTORIES
# ======================================

PARENT_DIR = "/home/marchandlab/DataAnalysis/Kaplan/basecall/8letter/251226_PB_Model-Testing/Mod"   # <-- EDIT THIS
MAX_ALIGNMENTS = 7500
RANDOM_SEED = 42

MASTER_RESULTS = []   # will store overall results from all child datasets


def define_directories(base_path):
    results_dir = os.path.join(base_path, f"alignment_results_capped_{MAX_ALIGNMENTS}")
    bed_dir = os.path.join(base_path, 'references')
    os.makedirs(results_dir, exist_ok=True)
    return results_dir, bed_dir

def read_bed_file(bed_dir):
    # Find any .bed file in the references directory
    bed_files = [f for f in os.listdir(bed_dir) if f.endswith(".bed")]

    if len(bed_files) == 0:
        print(f"No BED file found in: {bed_dir}")
        sys.exit(1)
    if len(bed_files) > 1:
        print(f"ERROR: Multiple BED files found in {bed_dir}: {bed_files}")
        print("This script expects exactly one BED file.")
        sys.exit(1)

    bed_path = os.path.join(bed_dir, bed_files[0])
    print(f"[INFO] Using BED file: {bed_path}")

    df = pd.read_csv(bed_path, delim_whitespace=True, header=None)
    df.columns = ['Alignment', 'Start', 'End', 'Name', 'Score', 'Strand']
    return df

def analyze_all_datasets(parent_dir):
    print(f"[INFO] Scanning parent directory: {parent_dir}")

    for entry in sorted(os.listdir(parent_dir)):
        base_path = os.path.join(parent_dir, entry)
        if not os.path.isdir(base_path):
            continue

        print(f"\n[INFO] Checking dataset: {base_path}")
        model_name = entry  # label for master results

        # Required folders/files
        references_dir = os.path.join(base_path, "references")
        remora_file = os.path.join(base_path, "remora_outputs/per-read_modifications.tsv")

        if not os.path.exists(references_dir):
            print(f"[SKIP] No references/ directory.")
            continue

        if not os.path.exists(remora_file):
            print(f"[SKIP] No per-read_modifications.tsv found.")
            continue

        # Process this dataset
        try:
            analyze_results(base_path)

            # After analysis → read the dataset's overall_results.csv
            overall_path = os.path.join(
                base_path,
                f"alignment_results_capped_{MAX_ALIGNMENTS}",
                "overall_results.csv"
            )

            if os.path.exists(overall_path):
                df_overall = pd.read_csv(overall_path)
                df_overall.insert(0, "Model", model_name)  # prepend dataset name
                MASTER_RESULTS.append(df_overall)
            else:
                print(f"[WARN] overall_results.csv missing for {model_name}")

        except Exception as e:
            print(f"[ERROR] Failed on {base_path}: {e}")

    # ---- After looping over all datasets: compile master CSV ----
    if MASTER_RESULTS:
        master_df = pd.concat(MASTER_RESULTS, ignore_index=True)
        out_path = os.path.join(parent_dir, "master_overall_results.csv")
        master_df.to_csv(out_path, index=False)
        print(f"\n[INFO] Master results written to: {out_path}")
    else:
        print("\n[INFO] No datasets produced overall results. Nothing to compile.")



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
        # ---- Check if below read cap ----
    if MAX_ALIGNMENTS is not None and len(df) < MAX_ALIGNMENTS:
        print(f"[WARN] Alignment {alignment} has only {len(df)} reads (< {MAX_ALIGNMENTS} cap)")

    df[['class_0_probs', 'class_1_probs']] = (
        df['class_probs'].str.split(',', expand=True).astype(float)
    )
    df['class_pred'] = (
        pd.to_numeric(df['class_pred'], errors='coerce')
        .fillna(0)
        .astype(int)
    )

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
    # ---- Apply capping here ----
    if MAX_ALIGNMENTS is not None:
        if len(df) < MAX_ALIGNMENTS:
            print(f"[WARN] Alignment {alignment} has only {len(df)} reads (< {MAX_ALIGNMENTS} cap)")
        elif len(df) > MAX_ALIGNMENTS:
            print(f"[INFO] Capping {alignment} from {len(df)} to {MAX_ALIGNMENTS} reads")
            df = df.sample(n=MAX_ALIGNMENTS, random_state=RANDOM_SEED)

        # Save individual alignment
        df.to_csv(full_alignment_results_path, mode='a',
                  header=not os.path.exists(full_alignment_results_path),
                  index=False)

        save_alignment_df(df, alignment, results_dir)
        result = calculate_results(df, alignment)
        full_results.append(result)



    # Save summary
    summary_df = pd.DataFrame(full_results)
    summary_df.to_csv(os.path.join(results_dir, 'overall_results.csv'), index=False)
    print(f"Analysis complete. Results saved to alignment_results_capped_{MAX_ALIGNMENTS}")

if __name__ == "__main__":
    analyze_all_datasets(PARENT_DIR)



