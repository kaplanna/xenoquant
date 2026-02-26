#!/usr/bin/env python3
# xr_signal_plot_step.py
# ------------------------------------------------------------
# Per-base step plot (median or trimmed mean) in reference coordinates.
# Compares canonical vs modified signal across aligned reads.
# ------------------------------------------------------------

from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pod5
import pysam
from remora import io, refine_signal_map
from xr_params import COL_STD, COL_MOD, COL_X_HIGHLIGHT, FLANK

# ============================================================
# === USER PARAMETERS ========================================
# ============================================================
WORKDIR_DEFAULT = Path("/path")
XNA_DEFAULT     = "S"
LEVELS_TXT      = "/path/to/9mer_10-4-1.tsv"
CHROM           = "contig1"
N_READS         = 500
ST_TRIM         = 1
EN_TRIM         = 1
LINE_W          = 2.0
ALPHA           = 0.95

# ============================================================
# === PLOTTING AESTHETICS ===================================
# ============================================================
mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 18,
    "axes.labelsize": 20,
    "axes.titlesize": 20,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 16,
    "axes.linewidth": 2.0,
    "lines.linewidth": 2.0,
    "xtick.major.width": 2.0,
    "ytick.major.width": 2.0,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "legend.frameon": False,
    "savefig.dpi": 300,
    "figure.dpi": 150,
    "axes.edgecolor": "black",
    "axes.grid": False,
})

# ============================================================
# === HELPER FUNCTIONS =======================================
# ============================================================

def style_axes(ax, xlabel=None, ylabel=None):
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)


def set_base_center_ticks(ax, start, end):
    """Auto-space x-axis tick labels based on region width."""
    centers = np.arange(start, end) + 0.5
    width = end - start
    if width <= 35:
        step = 1
    elif width <= 70:
        step = 2
    elif width <= 150:
        step = 5
    elif width <= 300:
        step = 10
    else:
        step = 20
    labels = [str(pos) if (i % step) == 0 else "" for i, pos in enumerate(range(start, end))]
    ax.set_xticks(centers)
    ax.set_xticklabels(labels, rotation=0, ha="center")


def set_standard_signal_yticks(ax):
    ax.set_yticks([-2, -1, 0, 1, 2])
    ax.set_yticklabels(["−2", "−1", "0", "+1", "+2"])



def load_xna_sites(bed_file, chrom):
    sites = []
    with open(bed_file) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.strip().split()
            if fields[0] == chrom:
                start, end = int(fields[1]), int(fields[2])
                name = fields[3] if len(fields) > 3 else "XNA"
                sites.append((start, end, name))
    return sites


def per_read_metric(sm_region, ref_start, ref_end, strand, metric="median", st_trim=1, en_trim=1):
    """Compute per-base signal metric (median or trimmed mean) for one read."""
    sig = sm_region.norm_signal
    s2s = sm_region.seq_to_sig_map
    nb = len(sm_region.seq)
    out = np.full(ref_end - ref_start, np.nan, float)

    def base_pos(b):
        return ref_start + b if strand == "+" else (ref_end - 1) - b

    for b in range(nb):
        s0, s1 = s2s[b], s2s[b + 1] if (b + 1) < len(s2s) else None
        if s0 is None or s1 is None or s1 <= s0:
            continue
        seg = sig[s0:s1]
        if seg.size == 0:
            continue
        if metric == "median":
            val = np.median(seg)
        elif metric == "trimmean":
            if seg.size <= (st_trim + en_trim):
                val = np.mean(seg)
            else:
                val = np.mean(seg[st_trim: seg.size - en_trim])
        else:
            raise ValueError("metric must be 'median' or 'trimmean'")
        idx = base_pos(b) - ref_start
        if 0 <= idx < len(out):
            out[idx] = val
    return out


def collect_metric_arrays(pod5_path, bam_path, chrom, start, end, n_reads, refiner,
                          metric="median", st_trim=1, en_trim=1):
    """Collect per-read signal metrics for all reads overlapping [start, end)."""
    dr = pod5.DatasetReader(pod5_path)
    bam = pysam.AlignmentFile(bam_path, "rb")

    per_read = []
    kept = total = 0
    for aln in bam.fetch(chrom, start, end):
        total += 1
        try:
            pod5_read = dr.get_read(aln.query_name)
            io_read = io.Read.from_pod5_and_alignment(pod5_read, aln)
            io_read.set_refine_signal_mapping(refiner, ref_mapping=False)
            io_read.set_refine_signal_mapping(refiner, ref_mapping=True)
        except Exception:
            continue

        sm_region = io_read.extract_ref_reg(
            io_read.ref_reg.adjust(
                start_adjust=start - io_read.ref_reg.start,
                end_adjust=end - io_read.ref_reg.end,
            )
        )
        arr = per_read_metric(sm_region, start, end, io_read.ref_reg.strand,
                              metric=metric, st_trim=st_trim, en_trim=en_trim)
        if np.any(np.isfinite(arr)):
            per_read.append(arr)
            kept += 1
        if kept >= n_reads:
            break

    arr = np.vstack(per_read) if kept else np.empty((0, end - start), float)
    print(f"[INFO] {Path(bam_path).stem}: kept={kept}/{total} reads.")
    return arr


# ============================================================
# === MAIN FUNCTION ==========================================
# ============================================================

def main():
    import argparse
    ap = argparse.ArgumentParser(description="Per-base step plot of normalized signal around an XNA site.")
    ap.add_argument("workdir", nargs="?", default=str(WORKDIR_DEFAULT))
    ap.add_argument("--xna", "-x", default=XNA_DEFAULT)
    ap.add_argument("--st-trim", type=int, default=ST_TRIM)
    ap.add_argument("--en-trim", type=int, default=EN_TRIM)
    args = ap.parse_args()

    workdir = Path(args.workdir).resolve()
    xna = args.xna.strip().upper()
    if len(xna) != 1 or not xna.isalpha():
        raise ValueError(f"[XNA] Must be a single letter A–Z. Got: {xna!r}")

    # --- Input setup ---
    can_root = workdir / "canonical_preprocess"
    mod_root = workdir / "modified_preprocess"
    can_pod5 = next((can_root / "pod5").glob("*.pod5"))
    mod_pod5 = next((mod_root / "pod5").glob("*.pod5"))
    can_bam  = can_root / "bam" / "aligned.sorted.bam"
    mod_bam  = mod_root / "bam" / "aligned.sorted.bam"
    xna_bed  = workdir / "references" / f"{xna}.bed"
    out_dir  = workdir / "signal_plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- Region ---
    sites = load_xna_sites(str(xna_bed), CHROM)
    if not sites:
        raise ValueError(f"No XNA sites found for {CHROM} in {xna_bed}")
    xna_start, xna_end, _ = sites[0]
    start, end = xna_start - FLANK, xna_end + FLANK

    refiner = refine_signal_map.SigMapRefiner(
        kmer_model_filename=LEVELS_TXT,
        do_rough_rescale=True,
        scale_iters=0,
        do_fix_guage=True,
    )

    # --- Compute both metrics ---
    results = {}
    for metric in ("median", "trimmean"):
        can_arr = collect_metric_arrays(str(can_pod5), str(can_bam), CHROM, start, end, N_READS, refiner,
                                        metric=metric, st_trim=args.st_trim, en_trim=args.en_trim)
        mod_arr = collect_metric_arrays(str(mod_pod5), str(mod_bam), CHROM, start, end, N_READS, refiner,
                                        metric=metric, st_trim=args.st_trim, en_trim=args.en_trim)
        results[metric] = (np.nanmean(can_arr, axis=0), np.nanmean(mod_arr, axis=0))

    # =========================================================
    # === PLOTTING ============================================
    # =========================================================
    for metric, (can_mean, mod_mean) in results.items():
        fig, ax = plt.subplots(figsize=(4.5, 3.25))
        positions = np.arange(start, end)
        # --- FIX: ensure last signal drawn ---
        edges = np.concatenate([positions, [end]])
        can_vals = np.concatenate([can_mean, [can_mean[-1]]])
        mod_vals = np.concatenate([mod_mean, [mod_mean[-1]]])

        ax.step(edges, can_vals, where="post", color=COL_STD, lw=LINE_W, alpha=ALPHA, label="Standard")
        ax.step(edges, mod_vals, where="post", color=COL_MOD, lw=LINE_W, alpha=ALPHA, label="Modified")
        ax.axvspan(xna_start, xna_start + 1, color=COL_X_HIGHLIGHT, alpha=0.2, lw=0)

        ax.set_xlim(start, end)
        set_standard_signal_yticks(ax)
        set_base_center_ticks(ax, start, end)
        style_axes(ax, xlabel=f"Position on {CHROM}", ylabel="Normalized signal")

        fig.tight_layout()
        out_stem = f"{workdir.name}_{xna}_{CHROM}_{start}-{end}_step_{metric}"
        for ext in ("svg", "pdf"):
            plt.savefig(out_dir / f"{out_stem}.{ext}")
        plt.close(fig)
        print(f"[INFO] Saved step plot: {out_stem}.svg/pdf")


if __name__ == "__main__":
    main()

