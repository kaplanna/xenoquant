
#!/usr/bin/env python3

import os
import subprocess
import sys

# ============================
# USER SETTINGS (EDIT THESE)
# ============================
FASTQ_DIR = "/home/marchandlab/DataAnalysis/Kaplan/training/8L/251203_8LT_training_demux"
OUTDIR_BASE = "/home/marchandlab/DataAnalysis/Kaplan/training/8L/251203_8LT_training_demux/nanoplot"
THREADS = 8

FASTQ_EXTENSIONS = (".fastq", ".fq", ".fastq.gz", ".fq.gz")

# ============================
# Sanity checks
# ============================
if not os.path.isdir(FASTQ_DIR):
    sys.exit(f"[ERROR] FASTQ directory not found: {FASTQ_DIR}")

os.makedirs(OUTDIR_BASE, exist_ok=True)

# ============================
# Find FASTQ files
# ============================
fastq_files = sorted(
    f for f in os.listdir(FASTQ_DIR)
    if f.endswith(FASTQ_EXTENSIONS)
)

if not fastq_files:
    sys.exit("[ERROR] No FASTQ files found")

print(f"[INFO] Found {len(fastq_files)} FASTQ files")

# ============================
# Run NanoPlot
# ============================
for fastq in fastq_files:
    fastq_path = os.path.join(FASTQ_DIR, fastq)

    # Clean sample name
    sample_name = fastq
    for ext in FASTQ_EXTENSIONS:
        if sample_name.endswith(ext):
            sample_name = sample_name[: -len(ext)]
            break

    outdir = os.path.join(OUTDIR_BASE, sample_name)
    os.makedirs(outdir, exist_ok=True)

    cmd = [
        "NanoPlot",
        "--fastq", fastq_path,
        "--outdir", outdir,
        "--threads", str(THREADS),
        "--loglength",
        "--N50"
    ]

    print(f"[INFO] Running NanoPlot for {fastq}")
    print("       " + " ".join(cmd))

    result = subprocess.run(cmd)

    if result.returncode != 0:
        print(f"[WARNING] NanoPlot failed for {fastq}")
    else:
        print(f"[INFO] Finished {fastq}")

print("[INFO] All NanoPlot runs completed")


