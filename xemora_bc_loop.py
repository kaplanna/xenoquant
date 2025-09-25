#!/usr/bin/env python3
"""
xemora_pipe_batch.py
Loop over a mixed_pod5 standard-curve directory, basecall each ratio/replicate,
then run results, demux, and raw basecall analysis per replicate (optional).

Tree expected (created by the mixer script):
MIXED_POD5_ROOT/
  fracB_0.00/rep1/*.pod5
  fracB_0.00/rep2/*.pod5
  ...
  fracB_1.00/rep3/*.pod5

Outputs will mirror this under BC_OUTPUT_ROOT.
"""

import os
import sys
import glob
import subprocess
from pathlib import Path
from typing import List

# -------------------------
# CONFIG — EDIT THESE
# -------------------------

# Root containing all ratios/replicates of mixed .pod5
MIXED_POD5_ROOT = Path("/home/marchandlab/DataAnalysis/Kaplan/basecall/Standard_Curve/mixed_pod5")

# Where to write basecall outputs (mirrors ratio/rep layout)
BC_OUTPUT_ROOT = Path("/home/marchandlab/DataAnalysis/Kaplan/basecall/Standard_Curve/ST_Curve")

# Reference/model for basecalling
BC_REF_FASTA = Path("/home/marchandlab/DataAnalysis/Kaplan/basecall/Standard_Curve/reference/ref_GBC.fasta")
BC_MODEL_FILE = Path("/home/marchandlab/github/kaplanna/xemora/models/240930_NTC_Models/GBC-ST-model_best.pt")

# Optional: demux mapping
BARCODE_PAIR_CSV = Path("/home/marchandlab/DataAnalysis/Kaplan/raw/xPCR/250814_P8_H4_dXTP_Opt/reference/DEMUX_P8_H4.csv")

# Which stages to run
RUN_BASECALL = True
RUN_ALIGNMENT_RESULTS = True
RUN_CUTADAPT_DEMUX = False
RUN_RAW_BASECALL_ANALYSIS = False

# Path to helper scripts (relative to this file, or absolute)
XR_RESULTS = Path("./lib/xr_results.py")
XR_DEMUX = Path("./lib/xr_demux.py")
XR_RAW_BC = Path("./lib/xr_raw_basecall_analysis.py")

# Python executable to use
PYTHON_EXE = sys.executable  # current interpreter
XEMORA_ENTRY = Path("xemora.py")  # adjust if it's elsewhere

# -------------------------
# IMPLEMENTATION
# -------------------------

def run(cmd: List[str]) -> None:
    """Run a command and stream output; raise on error."""
    print("[CMD]", " ".join(map(str, cmd)))
    subprocess.run(cmd, check=True)

def find_rep_dirs(root: Path) -> List[Path]:
    """
    Return all rep directories under MIXED_POD5_ROOT that contain at least one .pod5.
    e.g., [..., mixed_pod5/fracB_0.30/rep1, mixed_pod5/fracB_0.30/rep2, ...]
    """
    rep_dirs: List[Path] = []
    for ratio_dir in sorted([p for p in root.iterdir() if p.is_dir()]):
        for rep_dir in sorted([p for p in ratio_dir.iterdir() if p.is_dir()]):
            if any(rep_dir.glob("*.pod5")):
                rep_dirs.append(rep_dir)
            else:
                print(f"[WARN] No .pod5 found in {rep_dir} — skipping.")
    return rep_dirs

def main():
    if not MIXED_POD5_ROOT.exists():
        print(f"[ERROR] MIXED_POD5_ROOT not found: {MIXED_POD5_ROOT}")
        sys.exit(1)

    rep_dirs = find_rep_dirs(MIXED_POD5_ROOT)
    if not rep_dirs:
        print(f"[ERROR] No replicate directories with .pod5 files under {MIXED_POD5_ROOT}")
        sys.exit(2)

    print(f"[INFO] Found {len(rep_dirs)} replicate directories to process.")

    for rep_dir in rep_dirs:
        ratio_dir = rep_dir.parent.name         # e.g., 'fracB_0.30'
        rep_name = rep_dir.name                 # e.g., 'rep1'
        pod5_list = sorted(rep_dir.glob("*.pod5"))

        # Sanity
        if not pod5_list:
            print(f"[WARN] No .pod5 in {rep_dir}; skipping.")
            continue

        # Mirror output layout
        work_dir = BC_OUTPUT_ROOT / ratio_dir / rep_name
        work_dir.mkdir(parents=True, exist_ok=True)

        print("\n" + "="*90)
        print(f"[RUN ] Ratio={ratio_dir}  Rep={rep_name}")
        print(f"[INFO] Input dir: {rep_dir}")
        print(f"[INFO] Output (work_dir): {work_dir}")

        # ----------------- BASECALL -----------------
        if RUN_BASECALL:
            # Pass the replicate directory to -f so xemora finds the .pod5 inside
            cmd = [
                PYTHON_EXE, str(XEMORA_ENTRY), "basecall",
                "-w", str(work_dir),
                "-f", str(rep_dir),
                "-r", str(BC_REF_FASTA),
                "-m", str(BC_MODEL_FILE),
            ]
            run(cmd)

        # ----------------- ALIGNMENT RESULTS -----------------
        if RUN_ALIGNMENT_RESULTS:
            if not XR_RESULTS.exists():
                print(f"[WARN] xr_results.py not found at {XR_RESULTS}, skipping results.")
            else:
                cmd = [PYTHON_EXE, str(XR_RESULTS), str(work_dir)]
                run(cmd)

        # ----------------- DEMUX (optional) -----------------
        if RUN_CUTADAPT_DEMUX:
            if not XR_DEMUX.exists():
                print(f"[WARN] xr_demux.py not found at {XR_DEMUX}, skipping demux.")
            elif not BARCODE_PAIR_CSV.exists():
                print(f"[WARN] BARCODE_PAIR_CSV not found at {BARCODE_PAIR_CSV}, skipping demux.")
            else:
                cmd = [PYTHON_EXE, str(XR_DEMUX), str(work_dir), str(BARCODE_PAIR_CSV)]
                run(cmd)

        # ----------------- RAW BASECALL ANALYSIS (optional) -----------------
        if RUN_RAW_BASECALL_ANALYSIS:
            if not XR_RAW_BC.exists():
                print(f"[WARN] xr_raw_basecall_analysis.py not found at {XR_RAW_BC}, skipping raw analysis.")
            else:
                # FILTER_BY_CLASS_0 is usually imported from xr_params; pass as 'True'/'False' if needed.
                # Here we default to 'False'. Change to 'True' or pull from xr_params if desired.
                FILTER_BY_CLASS_0 = False
                cmd = [PYTHON_EXE, str(XR_RAW_BC), str(work_dir), str(int(FILTER_BY_CLASS_0))]
                run(cmd)

    print("\n[DONE] Batch basecalling finished for all ratio/replicate directories.")
    print(f"[OUTPUT ROOT] {BC_OUTPUT_ROOT}")

if __name__ == "__main__":
    main()

