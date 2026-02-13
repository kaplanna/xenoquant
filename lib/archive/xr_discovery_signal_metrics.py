#!/usr/bin/env python3
# xr_signal_metrics_discovery.py
# ----------------------------------------------------------------------
# Genome-wide per-base signal metrics comparison (canonical vs modified)
# No BED input required. Designed for exploratory modified-base detection
# on viral or small-genome datasets.
# ----------------------------------------------------------------------

from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pod5
import pysam
from Bio import SeqIO
from scipy.stats import ttest_ind
from remora import io, refine_signal_map
import csv

# ================================================================
# === USER CONFIGURATION =========================================
# ================================================================

# Canonical (PCR'd) dataset
CAN_POD5 = Path("/home/marchandlab/DataAnalysis/Kaplan/training/2509_Signal_Plots/250928_GBC_Plots/BA-plots/canonical_preprocess/pod5/10_pod5.pod5")
CAN_BAM  = Path("/home/marchandlab/DataAnalysis/Kaplan/training/2509_Signal_Plots/250928_GBC_Plots/BA-plots/canonical_preprocess/bam/aligned.sorted.bam")

# Modified (direct) dataset
MOD_POD5 = Path("/home/marchandlab/DataAnalysis/Kaplan/training/2509_Signal_Plots/250928_GBC_Plots/BA-plots/modified_preprocess/pod5/1_pod5.pod5")
MOD_BAM  = Path("/home/marchandlab/DataAnalysis/Kaplan/training/2509_Signal_Plots/250928_GBC_Plots/BA-plots/modified_preprocess/bam/aligned.sorted.bam")

# Reference FASTA (canonical sequence only)
REF_FASTA = Path("/home/marchandlab/DataAnalysis/Kaplan/training/2509_Signal_Plots/250928_GBC_Plots/BA-plots/references/xGBC_90mer_clean_N.fa")

# Model file for signal refinement (Remora 9-mer level table)
LEVELS_TXT = "/home/marchandlab/github/kaplanna/xemora/models/remora/9mer_10-4-1.tsv"

# Output directory
OUT_DIR = Path("/home/marchandlab/DataAnalysis/Kaplan/training/2509_Signal_Plots/250928_GBC_Plots/signal_discovery_testing")

# Number of reads to sample per condition
N_READS = 500

# Trim parameters for trimmed mean
ST_TRIM = 1
EN_TRIM = 1

# y-limits for stacked figure
YLIMS = {
    "trimmean": (-2, 2),
    "diff": (-1, 1),
    "dwell": (0, 40),
    "trimsd": (0, 0.4),
}

# ================================================================
# === GLOBAL STYLING =============================================
# ================================================================

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 11,
    "axes.labelsize": 11,
    "axes.titlesize": 11,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 11,
    "lines.linewidth": 2.0,
    "axes.linewidth": 2.0,
    "xtick.major.width": 2.0,
    "ytick.major.width": 2.0,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "savefig.dpi": 300,
    "figure.dpi": 150,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})

COL_STD = "#4575b4"         # canonical (blue)
COL_MOD = "#d73027"         # modified (red)

# ================================================================
# === FUNCTIONS ==================================================
# ================================================================

def style_axes(ax, xlabel=None, ylabel=None):
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)

def parse_reference(fasta_path):
    """Return (chrom, length) tuple for the first sequence in FASTA."""
    record = next(SeqIO.parse(str(fasta_path), "fasta"))
    return record.id, len(record.seq)

def refine_and_extract(dr, aln, refiner, start, end):
    """Return remora region + strand or None."""
    try:
        pod5_read = dr.get_read(aln.query_name)
        io_read = io.Read.from_pod5_and_alignment(pod5_read, aln)
        io_read.set_refine_signal_mapping(refiner, ref_mapping=False)
        io_read.set_refine_signal_mapping(refiner, ref_mapping=True)
        sm_region = io_read.extract_ref_reg(
            io_read.ref_reg.adjust(
                start_adjust=start - io_read.ref_reg.start,
                end_adjust=end - io_read.ref_reg.end))
        return sm_region, io_read.ref_reg.strand
    except Exception:
        return None, None

def per_base_metrics(sm_region, strand, region_len, st_trim=1, en_trim=1):
    """Return dict of per-base dwell, trimmean, trimsd."""
    sig = sm_region.norm_signal
    s2s = sm_region.seq_to_sig_map
    nb = len(sm_region.seq)
    dwell, trimmean, trimsd = [], [], []
    for b in range(nb):
        s0, s1 = s2s[b], s2s[b + 1]
        if s0 is None or s1 is None or s1 <= s0:
            continue
        seg = sig[s0:s1]
        if seg.size == 0:
            continue
        seg_t = seg[st_trim: seg.size - en_trim] if seg.size > (st_trim + en_trim) else seg
        dwell.append(seg.size)
        trimmean.append(np.mean(seg_t))
        trimsd.append(np.std(seg_t, ddof=1) if seg_t.size > 1 else 0.0)
    arrs = [np.array(a) for a in (dwell, trimmean, trimsd)]
    if strand == "-":
        arrs = [a[::-1] for a in arrs]
    return {"dwell": arrs[0], "trimmean": arrs[1], "trimsd": arrs[2]}

def collect_metrics(pod5_path, bam_path, chrom, ref_len, refiner, n_reads):
    """Collect per-base metrics across the entire reference."""
    dr = pod5.DatasetReader(pod5_path)
    bam = pysam.AlignmentFile(bam_path, "rb")
    dwell, trimmean, trimsd = [], [], []
    kept = total = 0
    
    
    for aln in bam.fetch(chrom):
        total += 1
        sm_reg, strand = refine_and_extract(dr, aln, refiner, 0, ref_len)
        if sm_reg is None:
            continue
        met = per_base_metrics(sm_reg, strand, ref_len)
        if met is None:
            continue
        for k in met:
            arr = met[k].astype(float)  # ensure float dtype for NaN padding
            met[k] = np.pad(arr, (0, max(0, ref_len - len(arr))),
                        mode='constant', constant_values=np.nan)[:ref_len]


        dwell.append(met["dwell"])
        trimmean.append(met["trimmean"])
        trimsd.append(met["trimsd"])
        kept += 1
        if kept >= n_reads:
            break

    mats = {
        "dwell": np.stack(dwell, axis=0) if dwell else np.empty((0, ref_len)),
        "trimmean": np.stack(trimmean, axis=0) if trimmean else np.empty((0, ref_len)),
        "trimsd": np.stack(trimsd, axis=0) if trimsd else np.empty((0, ref_len)),
    }
    print(f"[INFO] {Path(bam_path).stem}: kept={kept}/{total}")
    return mats

# ================================================================
# === MAIN =======================================================
# ================================================================

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    chrom, ref_len = parse_reference(REF_FASTA)
    print(f"[INFO] Reference: {chrom} (length={ref_len})")

    refiner = refine_signal_map.SigMapRefiner(
        kmer_model_filename=LEVELS_TXT, do_rough_rescale=True,
        scale_iters=0, do_fix_guage=True)

    can = collect_metrics(CAN_POD5, CAN_BAM, chrom, ref_len, refiner, N_READS)
    mod = collect_metrics(MOD_POD5, MOD_BAM, chrom, ref_len, refiner, N_READS)

    mean_can = {k: np.nanmean(can[k], axis=0) for k in can}
    mean_mod = {k: np.nanmean(mod[k], axis=0) for k in mod}
    diff_trim = mean_mod["trimmean"] - mean_can["trimmean"]

    # --- per-base Welch's t-test ---
    pvals = np.full(ref_len, np.nan)
    for i in range(ref_len):
        cvals = can["trimmean"][:, i]
        mvals = mod["trimmean"][:, i]
        if np.sum(np.isfinite(cvals)) >= 2 and np.sum(np.isfinite(mvals)) >= 2:
            pvals[i] = ttest_ind(cvals[np.isfinite(cvals)],
                                 mvals[np.isfinite(mvals)],
                                 equal_var=False).pvalue
    neglog10p = -np.log10(pvals)

    # --- write CSV summary ---
    csv_out = OUT_DIR / f"{chrom}_signal_metrics.csv"
    with open(csv_out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["position", "mean_can", "mean_mod", "diff_trim",
                    "dwell_can", "dwell_mod", "sd_can", "sd_mod", "p_value"])
        for i in range(ref_len):
            w.writerow([
                i + 1,
                mean_can["trimmean"][i],
                mean_mod["trimmean"][i],
                diff_trim[i],
                mean_can["dwell"][i],
                mean_mod["dwell"][i],
                mean_can["trimsd"][i],
                mean_mod["trimsd"][i],
                pvals[i],
            ])
    print(f"[INFO] Wrote per-base metrics to {csv_out}")

    # ============================================================
    # === PLOT ====================================================
    # ============================================================

    pos = np.arange(ref_len)
    fig, axes = plt.subplots(5, 1, figsize=(6.0, 14.0), sharex=True,
                             gridspec_kw={"hspace": 0.45})
    axA, axB, axC, axD, axE = axes

    def step(ax, y, color, lw=2.0):
        vals = np.concatenate([y, [y[-1]]])
        edges = np.concatenate([pos, [ref_len]])
        ax.step(edges, vals, where="post", lw=lw, color=color, alpha=0.95)

    # (A) Trimmed mean
    step(axA, mean_can["trimmean"], COL_STD)
    step(axA, mean_mod["trimmean"], COL_MOD)
    axA.set_ylabel("Trimmed mean")
    axA.set_ylim(*YLIMS["trimmean"])

    # (B) Δ Trimmed mean
    axB.axhline(0, color="#999999", linestyle="--", lw=1.5)
    step(axB, diff_trim, "#555555")
    axB.set_ylabel("Δ Trimmed mean")
    axB.set_ylim(*YLIMS["diff"])

    # (C) Dwell
    axC.scatter(pos + 0.5, mean_can["dwell"], color=COL_STD, s=16, alpha=0.9)
    axC.scatter(pos + 0.5, mean_mod["dwell"], color=COL_MOD, s=16, alpha=0.9)
    axC.set_ylabel("Dwell (samples)")
    axC.set_ylim(*YLIMS["dwell"])

    # (D) Trimmed SD
    axD.scatter(pos + 0.5, mean_can["trimsd"], color=COL_STD, s=16, alpha=0.9)
    axD.scatter(pos + 0.5, mean_mod["trimsd"], color=COL_MOD, s=16, alpha=0.9)
    axD.set_ylabel("Trimmed SD")
    axD.set_ylim(*YLIMS["trimsd"])

    # (E) −log10 p
    axE.scatter(pos + 0.5, neglog10p, color="#444444", s=20, alpha=0.8)
    axE.axhline(-np.log10(0.05), color="#A9A9A9", linestyle="--", lw=1.5)
    axE.set_ylabel("−log10 p")
    axE.set_xlabel("Reference position")

    for ax in axes:
        style_axes(ax)
        ax.margins(x=0)
        ax.set_xlim(0, ref_len)

    plt.subplots_adjust(top=0.97, bottom=0.08, left=0.12, right=0.97)
    pdf = OUT_DIR / f"{chrom}_signal_metrics.pdf"
    svg = OUT_DIR / f"{chrom}_signal_metrics.svg"
    plt.savefig(pdf)
    plt.savefig(svg)
    plt.show()
    plt.close(fig)
    print(f"[INFO] Wrote plots:\n - {pdf}\n - {svg}")

# ================================================================

if __name__ == "__main__":
    main()

