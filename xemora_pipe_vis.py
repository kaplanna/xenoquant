#!/usr/bin/env python3
########################################################################
# xemora_pipe.py
# ------------------------------------------------------------
# Master orchestration script for Xemora workflows:
#  - Model training
#  - Basecalling
#  - Alignment result extraction
#  - Cutadapt demultiplexing
#  - Raw basecall analysis
#  - Visualization suite (signal metrics, spaghetti, step, violin)
########################################################################

import os
import sys
import shlex
import subprocess
from pathlib import Path

from lib.xr_tools import *
from lib.xr_params import *

# ============================================================
# === CONFIGURATION ==========================================
# ============================================================

# --- Training paths ---
working_dir = '/home/marchandlab/DataAnalysis/Kaplan/training/2509_Signal_Plots/251226_GDsA_AT_Plots/DsA'
xna_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/DsPx/241209_G-Ds-A_Px-diol/20241209_1728_MN41475_AXC204_509514cb/1_pod5_test'
xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/DsPx/241209_G-Ds-A_Px-diol/reference/GDsA_60mer.fasta'
dna_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/DsPx/241210_GAA_CAT_xr_train/20241210_1413_MN37138_AWF246_a6bf4d1f/50_pod5_test'
dna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/DsPx/241209_G-Ds-A_Px-diol/reference/GDsA_60mer.fasta'

# --- Basecall paths ---
bc_working_dir = '/home/marchandlab/DataAnalysis/Kaplan/basecall/8letter/251219_H9-11_BC/BN_Ext-Train_5050_Basecall'
bc_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/xPCR/251218_H9-11_8L-Troubleshooting/20251218_1456_MN37138_BAN504_70d2149b/pod5'
bc_xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/xPCR/251218_H9-11_8L-Troubleshooting/reference/REF_H9-11_BN.fasta'
barcode_pair_csv = '/home/marchandlab/DataAnalysis/Kaplan/raw/xPCR/240816_B13_pH_Rerun_2/reference/DEMUX_B13.csv'
bc_model_file = '/home/marchandlab/DataAnalysis/Kaplan/training/8L/251207_8LT_ext_NN_train_demux/BN-Train/model/model_best.pt'

#bc_model_file = working_dir+'/model/model_best.pt'
# ============================================================
# === MASTER SWITCHES ========================================
# ============================================================

train_model              = False
basecall_reads           = True
output_alignment_results = False
cutadapt_demux           = False
raw_basecall_analysis    = True

# --- Visualization switches ---
plot_signal_metrics    = False   # calls xr_signal_metrics.py
plot_signal_spaghetti  = False   # calls xr_signal_plot_v2.py
plot_signal_step       = False   # calls xr_signal_plot_step.py
plot_signal_violin     = False    # calls xr_violin.py
extract_metrics        = False   # calls xr_extract_metrics.py

# ============================================================
# === HELPERS ================================================
# ============================================================

def get_xna_from_params() -> str:
    """Return XNA base (single uppercase letter) from xr_params."""
    xna = None
    try:
        xna = MOD_BASE  # type: ignore[name-defined]
    except Exception:
        pass
    if xna is None:
        try:
            xna = mod_base  # type: ignore[name-defined]
        except Exception:
            pass
    if not xna or not isinstance(xna, str) or len(xna.strip()) != 1 or not xna.strip().isalpha():
        raise ValueError("[XNA] Could not determine XNA (expected MOD_BASE or mod_base as single letter).")
    return xna.strip().upper()


def run(cmd_argv):
    """Print and run a subprocess command."""
    print("[RUN]", " ".join(shlex.quote(x) for x in cmd_argv))
    subprocess.run(cmd_argv, check=True)


# ============================================================
# === SCRIPT PATHS ===========================================
# ============================================================

HERE = Path(__file__).resolve().parent
PY = sys.executable or "python"

SCRIPT_METRICS   = str((HERE / "lib/xr_signal_metrics.py").resolve())
SCRIPT_SPAGHETTI = str((HERE / "lib/xr_signal_plot_v2.py").resolve())
SCRIPT_STEP      = str((HERE / "lib/xr_signal_plot_step.py").resolve())
SCRIPT_VIOLIN    = str((HERE / "lib/xr_violin.py").resolve())
SCRIPT_EXTRACT   = str((HERE / "lib/xr_extract_metrics.py").resolve())

# ============================================================
# === CORE PIPELINE ==========================================
# ============================================================

# --- 1. Train Xemora model ---
if train_model:
    cmd = [PY, "xemora.py", "train", "-w", working_dir,
           "-f", xna_fast5_dir, dna_fast5_dir,
           "-r", xna_ref_fasta, dna_ref_fasta]
    run(cmd)

# --- 2. Basecall reads ---
if basecall_reads:
    cmd = [PY, "xemora.py", "basecall",
           "-w", bc_working_dir,
           "-f", bc_fast5_dir,
           "-r", bc_xna_ref_fasta,
           "-m", bc_model_file]
    run(cmd)

# --- 3. Output alignment results ---
if output_alignment_results:
    cmd = [PY, str(HERE / "lib" / "xr_results.py"), bc_working_dir]
    run(cmd)

# --- 4. Cutadapt demultiplex ---
if cutadapt_demux:
    cmd = [PY, str(HERE / "lib" / "xr_demux.py"), bc_working_dir, barcode_pair_csv]
    run(cmd)

# --- 5. Raw basecall analysis ---
if raw_basecall_analysis:
    cmd = [PY, str(HERE / "lib" / "xr_raw_basecall_analysis.py"), bc_working_dir, str(FILTER_BY_CLASS_0)]
    run(cmd)

# ============================================================
# === 6. VISUALIZATIONS ======================================
# ============================================================

# If any visualization flag is True, run visualization stage
if any([plot_signal_metrics, plot_signal_spaghetti, plot_signal_step, plot_signal_violin, extract_metrics]):
    XNA = get_xna_from_params()
    WDIR = Path(working_dir).resolve()
    if not WDIR.exists():
        raise FileNotFoundError(f"[CONFIG] working_dir not found: {WDIR}")

    # Each visualization script accepts: <workdir> -x <XNA>
    if plot_signal_metrics:
        run([PY, SCRIPT_METRICS, str(WDIR), "-x", XNA])

    if plot_signal_spaghetti:
        run([PY, SCRIPT_SPAGHETTI, str(WDIR), "-x", XNA])

    if plot_signal_step:
        run([PY, SCRIPT_STEP, str(WDIR), "-x", XNA])

    if plot_signal_violin:
        run([PY, SCRIPT_VIOLIN, str(WDIR), "-x", XNA])

    if extract_metrics:
        run([PY, SCRIPT_EXTRACT, str(WDIR), "-x", XNA])

# ============================================================
# === END ====================================================
# ============================================================

