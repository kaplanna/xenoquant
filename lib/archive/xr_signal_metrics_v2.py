#!/usr/bin/env python3
# xr_signal_metrics_fullcov_pub.py
# Publication-ready per-base signal metrics with full-coverage enforcement.
# Includes optional ±SD ribbons (step-aligned) and boxplot dwell distributions.

import logging
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pod5
import polars as pl
import pysam
from scipy.stats import ttest_ind
from remora import io, refine_signal_map
from xr_params import COL_STD, COL_MOD, COL_X_HIGHLIGHT, FLANK

# ================================
# User parameters
# ================================
WORKDIR_DEFAULT = Path("/home/marchandlab/DataAnalysis/Kaplan/training/2509_Signal_Plots/250922_CNT_plots/Px-N")
XNA_DEFAULT     = "S"
CHROM = "contig1"
LEVELS_TXT = "/home/marchandlab/github/kaplanna/xemora/models/remora/9mer_10-4-1.tsv"
N_READS = 500

# ================================
# Display options
# ================================
SHOW_RIBBON = False  # Toggle ribbons globally

# Global y-limits
YLIMS = {
    "trimmean": (-2, 2),
    "diff": (-1, 1),
    "trimsd": (0, 0.4),
}

# ================================
# Global aesthetics
# ================================
mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 12,
    "axes.labelsize": 13,
    "axes.titlesize": 13,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
    "figure.titlesize": 14,
    "lines.linewidth": 2.0,
    "axes.linewidth": 2.0,
    "grid.linewidth": 1.0,
    "xtick.major.width": 2.0,
    "ytick.major.width": 2.0,
    "xtick.minor.width": 1.2,
    "ytick.minor.width": 1.2,
    "xtick.major.size": 4.5,
    "ytick.major.size": 4.5,
    "xtick.minor.size": 2.5,
    "ytick.minor.size": 2.5,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "legend.frameon": False,
    "savefig.dpi": 300,
    "figure.dpi": 150,
    "axes.edgecolor": "black",
    "axes.grid": False,
    "axes.axisbelow": False,
})

logging.getLogger("Remora").setLevel(logging.INFO)

# ---------- helpers ----------
def load_xna_sites(bed_file, chrom):
    with open(bed_file) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            flds = line.strip().split()
            if flds[0] == chrom:
                start, end = int(flds[1]), int(flds[2])
                name = flds[3] if len(flds) > 3 else "XNA"
                return start, end, name
    raise ValueError(f"No XNA sites found for {chrom} in {bed_file}")

def style_axes(ax, xlabel=None, ylabel=None):
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)

def refine_and_extract(dr, aln, refiner, start, end):
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

def fullcov_metrics(sm_region, strand, expected_len, st_trim=1, en_trim=1):
    sig = sm_region.norm_signal
    s2s = sm_region.seq_to_sig_map
    nb = len(sm_region.seq)
    if nb != expected_len or len(s2s) != nb + 1:
        return None
    dwell, mean, trimmn, trimsd = [], [], [], []
    for b in range(nb):
        s0, s1 = s2s[b], s2s[b + 1]
        if s0 is None or s1 is None or s1 <= s0:
            return None
        seg = sig[s0:s1]
        if seg.size == 0:
            return None
        seg_t = seg[st_trim: seg.size - en_trim] if seg.size > (st_trim + en_trim) else seg
        dwell.append(seg.size)
        mean.append(np.mean(seg))
        trimmn.append(np.mean(seg_t))
        trimsd.append(np.std(seg_t, ddof=1) if seg_t.size > 1 else 0.0)
    if strand == "-":
        dwell, mean, trimmn, trimsd = [np.array(arr[::-1]) for arr in (dwell, mean, trimmn, trimsd)]
    else:
        dwell, mean, trimmn, trimsd = [np.array(arr) for arr in (dwell, mean, trimmn, trimsd)]
    return {"dwell": dwell, "mean": mean, "trimmean": trimmn, "trimsd": trimsd}

def collect_sample_metrics(pod5_path, bam_path, chrom, start, end, n_reads, refiner):
    dr = pod5.DatasetReader(pod5_path)
    bam = pysam.AlignmentFile(bam_path, "rb")
    L = end - start
    mats = {"dwell": [], "mean": [], "trimmean": [], "trimsd": []}
    kept = total = 0
    for aln in bam.fetch(chrom, start, end):
        total += 1
        sm_reg, strand = refine_and_extract(dr, aln, refiner, start, end)
        if sm_reg is None:
            continue
        met = fullcov_metrics(sm_reg, strand, expected_len=L)
        if met is None:
            continue
        for k in mats:
            mats[k].append(met[k])
        kept += 1
        if kept >= n_reads:
            break
    for k in mats:
        mats[k] = np.stack(mats[k], axis=0) if mats[k] else np.empty((0, L), dtype=float)
    print(f"[INFO] {Path(bam_path).stem}: kept={kept}/{total} full coverage reads.")
    return mats

# ---------- main ----------
def main():
    import argparse
    ap = argparse.ArgumentParser(description="Publication-ready per-base signal metrics (full coverage only).")
    ap.add_argument("workdir", nargs="?", default=str(WORKDIR_DEFAULT))
    ap.add_argument("--xna", "-x", default=XNA_DEFAULT)
    args = ap.parse_args()

    workdir = Path(args.workdir).resolve()
    xna = args.xna.upper()

    can_root = workdir / "canonical_preprocess"
    mod_root = workdir / "modified_preprocess"
    can_pod5 = next((can_root / "pod5").glob("*.pod5"))
    mod_pod5 = next((mod_root / "pod5").glob("*.pod5"))
    can_bam  = can_root / "bam" / "aligned.sorted.bam"
    mod_bam  = mod_root / "bam" / "aligned.sorted.bam"
    xna_bed  = workdir / "references" / f"{xna}.bed"
    out_dir  = workdir / "signal_plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    xna_start, xna_end, _ = load_xna_sites(xna_bed, CHROM)
    start, end = xna_start - FLANK, xna_end + FLANK
    L = end - start
    pos = np.arange(start, end)

    refiner = refine_signal_map.SigMapRefiner(
        kmer_model_filename=LEVELS_TXT, do_rough_rescale=True,
        scale_iters=0, do_fix_guage=True)

    can = collect_sample_metrics(str(can_pod5), str(can_bam), CHROM, start, end, N_READS, refiner)
    mod = collect_sample_metrics(str(mod_pod5), str(mod_bam), CHROM, start, end, N_READS, refiner)

    # Check coverage
    cov_can = np.sum(np.isfinite(can["trimmean"]), axis=0)
    cov_mod = np.sum(np.isfinite(mod["trimmean"]), axis=0)
    if not np.all(cov_can == np.max(cov_can)) or not np.all(cov_mod == np.max(cov_mod)):
        raise RuntimeError("[ABORT] Incomplete coverage detected — stop before plotting.")

    mean_can = {k: np.nanmean(can[k], axis=0) for k in can}
    mean_mod = {k: np.nanmean(mod[k], axis=0) for k in mod}
    diff_trim = mean_mod["trimmean"] - mean_can["trimmean"]

    # Welch p-values
    pvals = np.array([
        ttest_ind(can["trimmean"][:, i][np.isfinite(can["trimmean"][:, i])],
                  mod["trimmean"][:, i][np.isfinite(mod["trimmean"][:, i])],
                  equal_var=False).pvalue
        if (np.sum(np.isfinite(can["trimmean"][:, i])) >= 2 and
            np.sum(np.isfinite(mod["trimmean"][:, i])) >= 2)
        else np.nan for i in range(L)
    ])
    with np.errstate(divide="ignore"):
        neglog10p = -np.log10(pvals)

    # ===============================
    # Single violin plot for trimmed mean signal with multi-level significance
    # ===============================

    P_THRESHOLDS = [(0.001, "***"), (0.01, "**"), (0.05, "*")]
    YLIM_TRIMMEAN = (-2, 2)

    def violin_plot_signal(metric_can, metric_mod, pvals, ylabel, ylim, out_suffix):
        fig, ax = plt.subplots(figsize=(6, 3.5))
        positions = np.arange(start, end)
        width = 0.35
        ymax_all = ylim[1]
        sig_y_gap = 0.05 * (ylim[1] - ylim[0])

        for i, p in enumerate(positions):
            can_vals = metric_can[:, i][np.isfinite(metric_can[:, i])]
            mod_vals = metric_mod[:, i][np.isfinite(metric_mod[:, i])]
            if len(can_vals) == 0 or len(mod_vals) == 0:
                continue

            # --- draw violins ---
            vdata = [can_vals, mod_vals]
            vpos = [p - width / 2, p + width / 2]
            parts = ax.violinplot(
                vdata,
                positions=vpos,
                widths=0.6,
                showmeans=False,
                showmedians=True,
                showextrema=False,
                bw_method=0.3
            )

            # Color violins
            for j, pc in enumerate(parts['bodies']):
                pc.set_facecolor([COL_STD, COL_MOD][j])
                pc.set_edgecolor([COL_STD, COL_MOD][j])
                pc.set_alpha(0.9)

            # Style median
            if 'cmedians' in parts:
                med = parts['cmedians']
                med.set_color('black')
                med.set_linewidth(1.2)

            # --- significance symbol ---
            pval = pvals[i]
            if np.isfinite(pval):
                sym = None
                for thr, mark in P_THRESHOLDS:
                    if pval < thr:
                        sym = mark
                        break
                if sym:
                    top_y = max(np.nanmax(can_vals), np.nanmax(mod_vals))
                    star_y = top_y + sig_y_gap
                    if star_y > ymax_all:
                        ymax_all = star_y + sig_y_gap
                    ax.text(p, star_y, sym, ha="center", va="bottom",
                            fontsize=14, fontweight="bold", color="black")

        # --- styling ---
        ax.axvspan(xna_start, xna_start + 1, color=COL_X_HIGHLIGHT, alpha=0.1)
        ax.set_xlim(start - 0.5, end + 0.5)
        ax.set_ylim(ylim[0], ymax_all)
        ax.set_xticks(positions)
        ax.set_xticklabels([str(p) for p in positions])
        ax.tick_params(axis="x", rotation=30)
        ax.set_ylabel(ylabel)
        style_axes(ax)
        plt.tight_layout()

        # --- save ---
        out_stem = f"{workdir.name}_{xna}_{CHROM}_{start}-{end}_{out_suffix}"
        pdf = out_dir / f"{out_stem}.pdf"
        svg = out_dir / f"{out_stem}.svg"
        plt.savefig(pdf, bbox_inches="tight")
        plt.savefig(svg, bbox_inches="tight")
        plt.show()
        plt.close(fig)
        print(f"[INFO] Saved violin plot: {pdf}")

    # === Generate ===
    violin_plot_signal(
        can["trimmean"], mod["trimmean"], pvals,
        ylabel="Trimmed mean signal", ylim=YLIM_TRIMMEAN,
        out_suffix="trimmean_violin"
    )



if __name__ == "__main__":
    main()

