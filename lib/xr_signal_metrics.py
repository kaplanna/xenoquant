#!/usr/bin/env python3
# xr_signal_metrics.py
# ------------------------------------------------------------
# Publication-ready per-base signal metrics with full coverage.
# Only + strand reads are analyzed (rev-comp handled upstream).
# Independent flank length (METRICS_FLANK); no ±SD ribbons.
# ------------------------------------------------------------

from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pod5
import pysam
from scipy.stats import ttest_ind
from remora import io, refine_signal_map
from xr_params import COL_STD, COL_MOD, COL_X_HIGHLIGHT

# ============================================================
# === USER PARAMETERS ========================================
# ============================================================
WORKDIR_DEFAULT = Path("/path")
XNA_DEFAULT     = "S"
CHROM           = "contig1"
LEVELS_TXT      = "/path/to/9mer_10-4-1.tsv"
N_READS         = 500
METRICS_FLANK   = 14  # independent flank length (bases) (21 for PZ & BS, 14 for DsPx due to short construct)

# Global y-limits
YLIMS = {
    "trimmean": (-2, 2),
    "diff": (-1, 1),
    "dwell": (0, 40),
    "trimsd": (0, 0.4),
}

# ============================================================
# === AESTHETICS =============================================
# ============================================================
mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 11,
    "axes.labelsize": 11,
    "axes.titlesize": 11,
    "xtick.labelsize": 11,
    "ytick.labelsize": 14,
    "legend.fontsize": 11,
    "figure.titlesize": 11,
    "lines.linewidth": 2.0,
    "axes.linewidth": 2.0,
    "xtick.major.width": 2,  
    "ytick.major.width": 2,   
    "xtick.direction": "out",
    "ytick.direction": "out",
    "legend.frameon": False,
    "savefig.dpi": 300,
    "figure.dpi": 150,
    "axes.edgecolor": "black",
    "axes.grid": False,
    "axes.axisbelow": False,
})




# ============================================================
# === HELPER FUNCTIONS =======================================
# ============================================================

def load_xna_sites(bed_file, chrom):
    """Read first matching site for chromosome."""
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

def set_base_center_ticks(ax, start, end):
    centers = np.arange(start, end) + 0.5
    width = end - start
    every = 1 if width <= 35 else 2 if width <= 70 else 5 if width <= 150 else 10
    labels = [str(p) if (i % every) == 0 else "" for i, p in enumerate(range(start, end))]
    ax.set_xticks(centers)
    ax.set_xticklabels(labels)

def refine_and_extract(dr, aln, refiner, start, end):
    """Refine signal mapping for one alignment and return region + strand."""
    try:
        pod5_read = dr.get_read(aln.query_name)
        io_read = io.Read.from_pod5_and_alignment(pod5_read, aln)
        io_read.set_refine_signal_mapping(refiner, ref_mapping=False)
        io_read.set_refine_signal_mapping(refiner, ref_mapping=True)
        sm_region = io_read.extract_ref_reg(
            io_read.ref_reg.adjust(
                start_adjust=start - io_read.ref_reg.start,
                end_adjust=end - io_read.ref_reg.end,
            )
        )
        return sm_region, io_read.ref_reg.strand
    except Exception:
        return None, None

def fullcov_metrics(sm_region, expected_len, st_trim=1, en_trim=1):
    """Compute dwell, mean, trimmed mean, trimmed SD for one + strand read."""
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
    return {
        "dwell": np.array(dwell),
        "mean": np.array(mean),
        "trimmean": np.array(trimmn),
        "trimsd": np.array(trimsd),
    }

def collect_sample_metrics(pod5_path, bam_path, chrom, start, end, n_reads, refiner):
    """Collect per-base metrics, excluding all - strand reads."""
    dr = pod5.DatasetReader(pod5_path)
    bam = pysam.AlignmentFile(bam_path, "rb")
    L = end - start
    mats = {k: [] for k in ("dwell", "mean", "trimmean", "trimsd")}
    kept = total = 0

    for aln in bam.fetch(chrom, start, end):
        total += 1
        sm_reg, strand = refine_and_extract(dr, aln, refiner, start, end)
        if sm_reg is None or strand == "-":
            continue  # skip reverse strand
        met = fullcov_metrics(sm_reg, expected_len=L)
        if met is None:
            continue
        for k in mats:
            mats[k].append(met[k])
        kept += 1
        if kept >= n_reads:
            break

    for k in mats:
        mats[k] = np.stack(mats[k], axis=0) if mats[k] else np.empty((0, L), float)
    print(f"[INFO] {Path(bam_path).stem}: kept={kept}/{total} (full coverage).")
    return mats

# ============================================================
# === MAIN ===================================================
# ============================================================

def main():
    import argparse
    ap = argparse.ArgumentParser(description="Per-base signal metrics (full coverage, + strand only).")
    ap.add_argument("workdir", nargs="?", default=str(WORKDIR_DEFAULT))
    ap.add_argument("--xna", "-x", default=XNA_DEFAULT)
    args = ap.parse_args()

    workdir = Path(args.workdir).resolve()
    xna = args.xna.upper()

    # Input paths
    can_root = workdir / "canonical_preprocess"
    mod_root = workdir / "modified_preprocess"
    can_pod5 = next((can_root / "pod5").glob("*.pod5"))
    mod_pod5 = next((mod_root / "pod5").glob("*.pod5"))
    can_bam  = can_root / "bam" / "aligned.sorted.bam"
    mod_bam  = mod_root / "bam" / "aligned.sorted.bam"
    xna_bed  = workdir / "references" / f"{xna}.bed"
    out_dir  = workdir / "signal_plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load target site
    xna_start, xna_end, _ = load_xna_sites(xna_bed, CHROM)
    start, end = xna_start - METRICS_FLANK, xna_end + METRICS_FLANK
    L = end - start
    pos = np.arange(start, end)
    centers = pos + 0.5

    refiner = refine_signal_map.SigMapRefiner(
        kmer_model_filename=LEVELS_TXT,
        do_rough_rescale=True,
        scale_iters=0,
        do_fix_guage=True,
    )

    # Collect metrics (+ strand only)
    can = collect_sample_metrics(str(can_pod5), str(can_bam), CHROM, start, end, N_READS, refiner)
    mod = collect_sample_metrics(str(mod_pod5), str(mod_bam), CHROM, start, end, N_READS, refiner)

    # Coverage check
    cov_can = np.sum(np.isfinite(can["trimmean"]), axis=0)
    cov_mod = np.sum(np.isfinite(mod["trimmean"]), axis=0)
    if not np.all(cov_can == np.max(cov_can)) or not np.all(cov_mod == np.max(cov_mod)):
        raise RuntimeError("[ABORT] Incomplete coverage detected — stop before plotting.")

    # Means + diffs
    mean_can = {k: np.nanmean(can[k], axis=0) for k in can}
    mean_mod = {k: np.nanmean(mod[k], axis=0) for k in mod}
    diff_trim = mean_mod["trimmean"] - mean_can["trimmean"]

    # p-values
    pvals = np.array([
        ttest_ind(
            can["trimmean"][:, i][np.isfinite(can["trimmean"][:, i])],
            mod["trimmean"][:, i][np.isfinite(mod["trimmean"][:, i])],
            equal_var=False,
        ).pvalue
        if (
            np.sum(np.isfinite(can["trimmean"][:, i])) >= 2
            and np.sum(np.isfinite(mod["trimmean"][:, i])) >= 2
        )
        else np.nan
        for i in range(L)
    ])
    with np.errstate(divide="ignore"):
        neglog10p = -np.log10(pvals)

    # =========================================================
    # === PLOTS ===============================================
    # =========================================================
    fig, axes = plt.subplots(
        nrows=5, ncols=1, figsize=(5, 14.0),
        gridspec_kw={"hspace": 0.45}, sharex=True
    )
    axA, axB, axC, axD, axE = axes
    edges = np.concatenate([pos, [end]])

    def step(ax, y, color):
        vals = np.concatenate([y, [y[-1]]])
        ax.step(edges, vals, where="post", lw=2.0, color=color, alpha=0.95)

    # (A) Trimmed mean
    axA.axvspan(xna_start, xna_start + 1, color=COL_X_HIGHLIGHT, alpha=0.2)
    step(axA, mean_can["trimmean"], COL_STD)
    step(axA, mean_mod["trimmean"], COL_MOD)
    axA.set_ylabel("Trimmed mean")
    axA.set_ylim(*YLIMS["trimmean"])

    # (B) Δ Trimmed mean
    axB.axvspan(xna_start, xna_start + 1, color=COL_X_HIGHLIGHT, alpha=0.2)
    axB.axhline(0.0, color="#A9A9A9", linestyle="--", lw=1.5)
    step(axB, diff_trim, "#555555")
    axB.set_ylabel("Δ Trimmed mean")
    axB.set_ylim(*YLIMS["diff"])

    # (C) Dwell
    axC.axvspan(xna_start, xna_start + 1, color=COL_X_HIGHLIGHT, alpha=0.2)
    axC.scatter(centers, mean_can["dwell"], color=COL_STD, s=18, alpha=0.9)
    axC.scatter(centers, mean_mod["dwell"], color=COL_MOD, s=18, alpha=0.9)
    axC.set_ylabel("Dwell (samples)")
    axC.set_ylim(*YLIMS["dwell"])

    # (D) Trimmed SD
    axD.axvspan(xna_start, xna_start + 1, color=COL_X_HIGHLIGHT, alpha=0.2)
    axD.scatter(centers, mean_can["trimsd"], color=COL_STD, s=18, alpha=0.9)
    axD.scatter(centers, mean_mod["trimsd"], color=COL_MOD, s=18, alpha=0.9)
    axD.set_ylabel("Trimmed SD")
    axD.set_ylim(*YLIMS["trimsd"])

    # (E) -log10 p-value
    axE.axvspan(xna_start, xna_start + 1, color=COL_X_HIGHLIGHT, alpha=0.2)
    axE.scatter(centers, neglog10p, color="#333333", s=22, alpha=0.8)
    axE.axhline(-np.log10(0.05), color="#A9A9A9", linestyle="--", lw=1.5, label="p = 0.05")
    axE.set_ylabel("−log10 p")
    axE.set_xlabel("Reference position")
    axE.legend(frameon=False, loc="upper right")

    # Style + save
    for ax in axes:
        ax.margins(x=0)
        ax.set_xlim(start, end)
        style_axes(ax)
    set_base_center_ticks(axE, start, end)
    plt.subplots_adjust(top=0.97, bottom=0.08, left=0.14, right=0.97)

    out_stem = f"{workdir.name}_{xna}_{CHROM}_{start}-{end}_metrics"
    for ext in ("pdf", "svg"):
        plt.savefig(out_dir / f"{out_stem}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Saved metrics plots: {out_stem}.pdf/svg")

if __name__ == "__main__":
    main()

