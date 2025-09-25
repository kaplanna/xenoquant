########################################################################
########################################################################
"""
xemora_pipe.py
"""
########################################################################
########################################################################

import os
import subprocess
import sys
import glob
import shlex
from pathlib import Path

from lib.xr_tools import *
from lib.xr_params import *  # expects MOD_BASE or mod_base to exist

############################################################
# Training paths
working_dir = '/home/marchandlab/DataAnalysis/Kaplan/training/2509_Signal_Plots/250925_PZ_71mer_Training_plot/PG-Train'
xna_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/PZ/250905_PZ_71mer_xr_train/20250905_1526_MN41475_BAC309_c46e04e9/pod5_skip'
xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/PZ/250905_PZ_71mer_xr_train/reference/71mer_PZ.fasta'
dna_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/PZ/240216_GC_71merPCR_xr_Train/20240216_1817_MN41475_ASE526_f9fc38c7/300_pod5_train'
dna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/PZ/250905_PZ_71mer_xr_train/reference/71mer_PZ.fasta'

############################################################
# Basecall paths
bc_working_dir = '/home/marchandlab/DataAnalysis/Kaplan/basecall/8letter/250921_P8_H4_NN_Basecall/SN-Basecall'
bc_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/xPCR/250814_P8_H4_dXTP_Opt/20250814_1746_MN37138_AYK921_bdc73555/pod5'
bc_xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/xPCR/250814_P8_H4_dXTP_Opt/reference/REF_P8_H4_NB.fasta'
barcode_pair_csv = '/home/marchandlab/DataAnalysis/Kaplan/raw/xPCR/250814_P8_H4_dXTP_Opt/reference/DEMUX_P8_H4.csv'
bc_model_file = '/home/marchandlab/github/kaplanna/xemora/models/250901_GNT_BSn_Models/SN-model_best.pt'
# bc_model_file = working_dir+'/model/model_best.pt'

############################################################
# Master switches
train_model = True
basecall_reads = False
output_alignment_results = False
cutadapt_demux = False
raw_basecall_analysis = False

# NEW: visualization switches (run after training outputs exist)
plot_signal_metrics   = True   # calls xr_signal_metrics.py
plot_signal_spaghetti = True   # calls xr_signal_plot_v2.py
plot_signal_step      = True   # calls xr_signal_plot_step.py

############################################################
# Helpers

def get_xna_from_params() -> str:
    """Prefer MOD_BASE (upper) then mod_base (lower). Must be a single letter."""
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
        raise ValueError("[XNA] Could not determine XNA from xr_params (expected MOD_BASE or mod_base as one letter).")
    return xna.strip().upper()

def run(cmd_argv):
    print("[RUN]", " ".join(shlex.quote(x) for x in cmd_argv))
    subprocess.run(cmd_argv, check=True)


# Resolve scripts (assume they live next to this file)
HERE = Path(__file__).resolve().parent
SCRIPT_METRICS   = str((HERE/ "lib/xr_signal_metrics.py").resolve())
SCRIPT_SPAGHETTI = str((HERE / "lib/xr_signal_plot_v2.py").resolve())
SCRIPT_STEP      = str((HERE / "lib/xr_signal_plot_step.py").resolve())

PY = sys.executable or "python"

############################################################
# Train dataset with Xemora train
if train_model:
    cmd = [PY, "xemora.py", "train", "-w", working_dir, "-f", xna_fast5_dir, dna_fast5_dir, "-r", xna_ref_fasta, dna_ref_fasta]
    run(cmd)

# Basecall fast5 directory 
if basecall_reads:
    cmd = [PY, "xemora.py", "basecall", "-w", bc_working_dir, "-f", bc_fast5_dir, "-r", bc_xna_ref_fasta, "-m", bc_model_file]
    run(cmd)

# Output results
if output_alignment_results:
    results_path = str(HERE / "lib" / "xr_results.py")
    cmd = [PY, results_path, bc_working_dir]
    run(cmd)

# Demux
if cutadapt_demux:
    demux_path = str(HERE / "lib" / "xr_demux.py")
    cmd = [PY, demux_path, bc_working_dir, barcode_pair_csv]
    run(cmd)

# Raw basecall analysis
if raw_basecall_analysis:
    raw_bc_path = str(HERE / "lib" / "xr_raw_basecall_analysis.py")
    cmd = [PY, raw_bc_path, bc_working_dir, str(FILTER_BY_CLASS_0)]
    run(cmd)

############################################################
# Visualizations (training WORKDIR + XNA from params)
if any([plot_signal_metrics, plot_signal_spaghetti, plot_signal_step]):
    XNA = get_xna_from_params()
    WDIR = Path(working_dir).resolve()
    if not WDIR.exists():
        raise FileNotFoundError(f"[CONFIG] working_dir not found: {WDIR}")

    # Each viz script now accepts: <workdir> -x <XNA> and writes to {workdir}/signal_plots/
    if plot_signal_metrics:
        # metrics (full-coverage, per-base stats and 5-panel plot)
        cmd = [PY, SCRIPT_METRICS, str(WDIR), "-x", XNA]
        run(cmd)

    if plot_signal_spaghetti:
        # spaghetti/polyline plot (compressed time by base bins)
        cmd = [PY, SCRIPT_SPAGHETTI, str(WDIR), "-x", XNA]
        run(cmd)

    if plot_signal_step:
        # per-base step plot (median/trimmean). You can add CLI options here if you want:
        # e.g., "--metric", "trimmean", "--st-trim", "1", "--en-trim", "1"
        cmd = [PY, SCRIPT_STEP, str(WDIR), "-x", XNA]
        run(cmd)

