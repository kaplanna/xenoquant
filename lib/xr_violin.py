#!/usr/bin/env python3
# xr_violin.py
# ------------------------------------------------------------
# Publication-ready violin plots for per-base trimmed-mean signals.
# ------------------------------------------------------------
# Workflow:
#   1. Load canonical & modified datasets (pod5 + BAM)
#   2. Refine signal mapping (Remora refiner)
#   3. Compute per-base trimmed mean signal (full coverage)
#   4. Perform Welch’s t-test per base
#   5. Generate violin plots with significance stars
# ------------------------------------------------------------

import logging
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pod5
import pysam
from scipy.stats import ttest_ind
from remora import io, refine_signal_map
from xr_params import COL_STD, COL_MOD, COL_X_HIGHLIGHT, FLANK

# ============================================================
# === USER PARAMETERS ========================================
# ============================================================
WORKDIR_DEFAULT = Path(
    ""
)
XNA_DEFAULT = "S"
CHROM = "contig1"
LEVELS_TXT = "/path/to/9mer_10-4-1.tsv"
N_READS = 500

# ============================================================
# === MATPLOTLIB AESTHETICS =================================
# ============================================================
mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 14,
        "axes.labelsize": 16,
        "axes.titlesize": 16,
        "xtick.labelsize": 12,
        "ytick.labelsize": 16,
        "legend.fontsize": 12,
        "figure.titlesize": 14,
        "lines.linewidth": 2.0,
        "axes.linewidth": 2.0,
        "grid.linewidth": 1.0,
        "xtick.major.width": 2.0,
        "ytick.major.width": 2.0,
        "legend.frameon": False,
        "savefig.dpi": 300,
        "figure.dpi": 150,
        "axes.edgecolor": "black",
        "axes.grid": False,
    }
)

logging.getLogger("Remora").setLevel(logging.INFO)

# ============================================================
# === HELPER FUNCTIONS =======================================
# ============================================================

def load_xna_sites(bed_file, chrom):
    """Return (start, end, name) for the XNA region from BED file."""
    with open(bed_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            flds = line.strip().split()
            if flds[0] == chrom:
                start, end = int(flds[1]), int(flds[2])
                name = flds[3] if len(flds) > 3 else "XNA"
                return start, end, name
    raise ValueError(f"No XNA sites found for {chrom} in {bed_file}")


def style_axes(ax, xlabel=None, ylabel=None):
    """Remove top/right spines and apply uniform labeling."""
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)


def refine_and_extract(dr, aln, refiner, start, end):
    """
    Refine signal mapping for a single read and extract normalized
    signal segment corresponding to [start, end).
    """
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


def fullcov_metrics(sm_region, strand, expected_len, st_trim=1, en_trim=1):
    """
    Compute per-base dwell, mean, trimmed mean, and trimmed SD.
    Returns dict of arrays if full coverage is valid, else None.
    """
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
        seg_t = seg[st_trim : seg.size - en_trim] if seg.size > (st_trim + en_trim) else seg
        dwell.append(seg.size)
        mean.append(np.mean(seg))
        trimmn.append(np.mean(seg_t))
        trimsd.append(np.std(seg_t, ddof=1) if seg_t.size > 1 else 0.0)

    arrays = [np.array(arr[::-1]) if strand == "-" else np.array(arr) for arr in (dwell, mean, trimmn, trimsd)]
    return dict(zip(["dwell", "mean", "trimmean", "trimsd"], arrays))


def collect_sample_metrics(pod5_path, bam_path, chrom, start, end, n_reads, refiner):
    """
    Iterate through BAM alignments overlapping [start, end).
    Collect per-base trimmed mean and other metrics across full-coverage reads.
    """
    dr = pod5.DatasetReader(pod5_path)
    bam = pysam.AlignmentFile(bam_path, "rb")
    L = end - start
    mats = {k: [] for k in ("dwell", "mean", "trimmean", "trimsd")}
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
        mats[k] = np.stack(mats[k], axis=0) if mats[k] else np.empty((0, L), float)
    print(f"[INFO] {Path(bam_path).stem}: kept={kept}/{total} full coverage reads.")
    return mats


# ============================================================
# === MAIN PIPELINE ==========================================
# ============================================================

def main():
    """Run full pipeline: collect metrics → compute stats → plot violins."""
    import argparse

    ap = argparse.ArgumentParser(description="Generate violin plots for per-base trimmed mean signal.")
    ap.add_argument("workdir", nargs="?", default=str(WORKDIR_DEFAULT))
    ap.add_argument("--xna", "-x", default=XNA_DEFAULT)
    args = ap.parse_args()

    workdir = Path(args.workdir).resolve()
    xna = args.xna.upper()

    # --- I/O paths ---
    can_root = workdir / "canonical_preprocess"
    mod_root = workdir / "modified_preprocess"
    can_pod5 = next((can_root / "pod5").glob("*.pod5"))
    mod_pod5 = next((mod_root / "pod5").glob("*.pod5"))
    can_bam = can_root / "bam" / "aligned.sorted.bam"
    mod_bam = mod_root / "bam" / "aligned.sorted.bam"
    xna_bed = workdir / "references" / f"{xna}.bed"
    out_dir = workdir / "signal_plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- Reference region setup ---
    xna_start, xna_end, _ = load_xna_sites(xna_bed, CHROM)
    start, end = xna_start - FLANK, xna_end + FLANK
    L = end - start

    # --- Initialize Remora refiner ---
    refiner = refine_signal_map.SigMapRefiner(
        kmer_model_filename=LEVELS_TXT,
        do_rough_rescale=True,
        scale_iters=0,
        do_fix_guage=True,
    )

    # --- Collect canonical / modified data ---
    can = collect_sample_metrics(str(can_pod5), str(can_bam), CHROM, start, end, N_READS, refiner)
    mod = collect_sample_metrics(str(mod_pod5), str(mod_bam), CHROM, start, end, N_READS, refiner)

    # --- Validate full coverage ---
    cov_can = np.sum(np.isfinite(can["trimmean"]), axis=0)
    cov_mod = np.sum(np.isfinite(mod["trimmean"]), axis=0)
    if not np.all(cov_can == np.max(cov_can)) or not np.all(cov_mod == np.max(cov_mod)):
        raise RuntimeError("[ABORT] Incomplete coverage detected — stop before plotting.")

    # --- Per-base Welch’s t-test ---
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

    # =========================================================
    # === VIOLIN PLOT =========================================
    # =========================================================
    P_THRESHOLDS = [(0.001, "***"), (0.01, "**"), (0.05, "*")]

    def violin_plot_signal(metric_can, metric_mod, pvals, ylabel, out_suffix):
        """Draw paired violins per base, annotate significance."""
        fig, ax = plt.subplots(figsize=(6, 3.5))
        positions = np.arange(start, end)
        width = 0.35
        sig_y_gap = 0.1
        ymax_all = 2

        for i, p in enumerate(positions):
            can_vals = metric_can[:, i][np.isfinite(metric_can[:, i])]
            mod_vals = metric_mod[:, i][np.isfinite(metric_mod[:, i])]
            if len(can_vals) == 0 or len(mod_vals) == 0:
                continue

            # --- Draw paired violins ---
            vdata = [can_vals, mod_vals]
            vpos = [p - width / 2, p + width / 2]
            parts = ax.violinplot(
                vdata,
                positions=vpos,
                widths=0.6,
                showmeans=False,
                showmedians=True,
                showextrema=False,
                bw_method=0.3,
            )

            # Color violins
            for j, pc in enumerate(parts["bodies"]):
                pc.set_facecolor([COL_STD, COL_MOD][j])
                pc.set_edgecolor([COL_STD, COL_MOD][j])
                pc.set_alpha(0.9)

            # Median line
            if "cmedians" in parts and parts["cmedians"] is not None:
                parts["cmedians"].set_color("white")
                parts["cmedians"].set_linewidth(1.2)

            # --- Add significance stars ---
            pval = pvals[i]
            if np.isfinite(pval):
                sym = next((mark for thr, mark in P_THRESHOLDS if pval < thr), None)
                if sym:
                    top_y = max(np.nanmax(can_vals), np.nanmax(mod_vals))
                    star_y = top_y + sig_y_gap
                    ymax_all = max(ymax_all, star_y + sig_y_gap)
                    ax.text(p, star_y, sym, ha="center", va="bottom",
                            fontsize=14, fontweight="bold", color="black")

        # --- Style & annotate XNA site ---
        ax.axvspan(xna_start - 0.5, xna_start + 0.5, color=COL_X_HIGHLIGHT, alpha=0.1)
        ax.set_xlim(start - 0.5, end + 0.5)
        ax.set_ylim(-3, 3)
        ax.set_yticks(np.arange(-3, 3.5, 3))
        ax.set_xticks(positions)
        ax.set_xticklabels([str(p) for p in positions], rotation=30)
        ax.set_ylabel(ylabel)
        style_axes(ax)
        plt.tight_layout()

        # --- Save outputs ---
        out_stem = f"{workdir.name}_{xna}_{CHROM}_{start}-{end}_{out_suffix}"
        for ext in ("pdf", "svg"):
            plt.savefig(out_dir / f"{out_stem}.{ext}", bbox_inches="tight")

        plt.close(fig)
        print(f"[INFO] Saved violin plot: {out_stem}.pdf/svg")

    # === Run plot ===
    violin_plot_signal(
        can["trimmean"],
        mod["trimmean"],
        pvals,
        ylabel="Trimmed mean signal",
        out_suffix="trimmean_violin",
    )


if __name__ == "__main__":
    main()

