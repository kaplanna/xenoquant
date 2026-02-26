#!/usr/bin/env python3
# xr_levels_spaghetti_lines_pub.py
# Remora-style "compressed-time" lines per read on base coordinates (strand-aware), publication styling.

from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
import pod5
import pysam
from remora import io, refine_signal_map
from xr_params import COL_STD, COL_MOD, COL_X_HIGHLIGHT, FLANK



LEVELS_TXT = "/path/to/9mer_10-4-1.tsv"
CHROM = "contig1"
N_READS = 75


# Rendering mode: "line" (compressed-time polyline), "dots" (scatter), "micro" (tiny segments)
RENDER_MODE = "line"


# Plot tuning
POINT_SIZE  = 2.0
ALPHA_STD   = 0.15
ALPHA_MOD   = 0.15
BIN_JITTER  = 0.08     # 0..0.5 spreads samples within base bin
LABEL_EVERY = None     # auto if None
SHOW_X_LABEL = True
SHOW_Y_LABEL = True
SHOW_XNA_CENTER_LINE = False

# ================================
# Global aesthetics
# ================================
mpl.rcParams.update({
    # ---- Fonts ----
    "font.family": "Arial",
    "font.size": 18,          # base font size
    "axes.labelsize": 20,     # axis labels
    "axes.titlesize": 20,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 16,
    "figure.titlesize": 22,

    # ---- Line and tick weights ----
    "axes.linewidth": 2.0,    # main spines
    "lines.linewidth": 2.0,   # general line plots
    "xtick.major.width": 2.0,
    "ytick.major.width": 2.0,
    "xtick.minor.width": 1.5,
    "ytick.minor.width": 1.5,
    "xtick.major.size": 6.0,
    "ytick.major.size": 6.0,
    "xtick.minor.size": 3.0,
    "ytick.minor.size": 3.0,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.pad": 6.0,
    "ytick.major.pad": 6.0,

    # ---- Figure and legend ----
    "legend.frameon": False,
    "legend.handlelength": 1.8,
    "legend.handletextpad": 0.5,
    "axes.labelpad": 10,
    "axes.titlepad": 12,

    # ---- Save and rendering ----
    "savefig.dpi": 300,
    "figure.dpi": 150,
    "savefig.transparent": False,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,

    # ---- Spines and grid ----
    "axes.edgecolor": "black",
    "axes.grid": False,
    "axes.axisbelow": False,
    "text.usetex": False,
})

def style_axes(ax, *, xlabel=None, ylabel=None):
    """Keep only left/bottom axes visible, no enclosing box, with y-label retained."""
    # Show only left and bottom spines (x and y axes)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_visible(True)
        ax.spines[spine].set_linewidth(2.0)

    # Keep default background; just avoid extra frame
    ax.set_frame_on(True)

    # Labels
    if SHOW_X_LABEL and xlabel:
        ax.set_xlabel(xlabel, labelpad=8)
    if SHOW_Y_LABEL and ylabel:
        ax.set_ylabel(ylabel, labelpad=8)

def set_standard_signal_yticks(ax):
    """Force y-axis ticks and labels to [-2, 0, +2] for signal-style plots."""
    ax.set_yticks([-2, 0, 2])
    ax.set_yticklabels(["−2", "0", "+2"])




def set_base_center_ticks(ax, start, end, *, label_every=None):
    centers = np.arange(start, end) + 0.5
    width = end - start
    if label_every is None:
        if width <= 35:      label_every = 1
        elif width <= 70:    label_every = 2
        elif width <= 150:   label_every = 5
        else:                label_every = 10
    labels = [str(p) if (k % label_every) == 0 else "" for k, p in enumerate(range(start, end))]
    ax.set_xticks(centers)
    ax.set_xticklabels(labels)

def load_xna_sites(bed_file, chrom):
    sites = []
    with open(bed_file) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.strip().split()
            if fields[0] == chrom:
                start = int(fields[1]); end = int(fields[2])
                name = fields[3] if len(fields) > 3 else "XNA"
                sites.append((start, end, name))
    return sites

# ---------- Render helpers ----------

def compressed_time_polyline(sm_region, ref_start, ref_end, strand, jitter=0.08):
    sig = sm_region.norm_signal
    s2s = sm_region.seq_to_sig_map
    nb  = len(sm_region.seq)
    n_sig = len(sig)

    def base_left(b):
        return (ref_start + b) if strand == "+" else ((ref_end - 1) - b)

    xs_parts, ys_parts = [], []
    left  = max(0.0, min(0.5, jitter))
    right = 1.0 - left

    for b in range(nb):
        s0 = s2s[b]
        s1 = s2s[b+1] if (b+1) < len(s2s) else None
        if s0 is None or s1 is None or s0 < 0 or s1 > n_sig or s1 <= s0:
            xs_parts.append(np.array([np.nan])); ys_parts.append(np.array([np.nan])); continue
        seg = sig[s0:s1]; n = seg.size
        L = base_left(b)
        xs = np.array([L + 0.5]) if n == 1 else (L + np.linspace(left, right, n))
        xs_parts.append(xs); ys_parts.append(seg.astype(float, copy=False))

    if not xs_parts:
        return np.empty(0), np.empty(0)
    return np.concatenate(xs_parts), np.concatenate(ys_parts)

def micro_segments(ax, xs, ys, color, alpha=0.15, seg_frac=0.6, lw=0.4):
    if xs.size == 0: return
    half = 0.5 * seg_frac
    x0 = xs - half; x1 = xs + half
    segs = np.stack([np.column_stack([x0, ys]), np.column_stack([x1, ys])], axis=1)
    lc = LineCollection(segs, colors=[color], linewidths=lw, alpha=alpha)
    ax.add_collection(lc)

def samples_as_points(sm_region, ref_start, ref_end, strand, jitter=0.10):
    sig = sm_region.norm_signal
    s2s = sm_region.seq_to_sig_map
    nb  = len(sm_region.seq)
    n_sig = len(sig)

    def base_left(b):
        return (ref_start + b) if strand == "+" else ((ref_end - 1) - b)

    left  = max(0.0, min(0.5, jitter))
    right = 1.0 - left
    xs_list, ys_list = [], []
    for b in range(nb):
        s0 = s2s[b]
        s1 = s2s[b+1] if (b+1) < len(s2s) else None
        if s0 is None or s1 is None or s0 < 0 or s1 > n_sig or s1 <= s0:
            continue
        seg = sig[s0:s1]; n = seg.size
        L = base_left(b)
        xs = np.array([L + 0.5]) if n == 1 else (L + np.linspace(left, right, n))
        xs_list.append(xs); ys_list.append(seg.astype(float, copy=False))
    if not xs_list: return np.empty(0), np.empty(0)
    return np.concatenate(xs_list), np.concatenate(ys_list)

# ---------- Collection ----------

def collect_lines(pod5_path, bam_path, chrom, start, end, n_reads, refiner):
    dr = pod5.DatasetReader(pod5_path)
    bam_fh = pysam.AlignmentFile(bam_path, "rb")
    lines = []; kept = total = 0
    for bam_aln in bam_fh.fetch(chrom, start, end):
        total += 1
        try:
            pod5_read = dr.get_read(bam_aln.query_name)
        except KeyError:
            continue

        io_read = io.Read.from_pod5_and_alignment(pod5_read, bam_aln)
        io_read.set_refine_signal_mapping(refiner, ref_mapping=False)
        io_read.set_refine_signal_mapping(refiner, ref_mapping=True)

        sm_region = io_read.extract_ref_reg(
            io_read.ref_reg.adjust(
                start_adjust=start - io_read.ref_reg.start,
                end_adjust=end   - io_read.ref_reg.end
            )
        )
        strand = io_read.ref_reg.strand

        if RENDER_MODE == "line":
            xs, ys = compressed_time_polyline(sm_region, start, end, strand, jitter=BIN_JITTER)
        elif RENDER_MODE in ("dots", "micro"):
            xs, ys = samples_as_points(sm_region, start, end, strand, jitter=BIN_JITTER)
        else:
            raise ValueError(f"Unknown RENDER_MODE: {RENDER_MODE}")

        if xs.size:
            lines.append((xs, ys))
            kept += 1
            if kept >= n_reads:
                break
    return lines, kept, total

# ---------- Main (all path logic here) ----------


def main():
    import argparse

    ap = argparse.ArgumentParser(
        description="Spaghetti lines of normalized signal around an XNA site."
    )
    ap.add_argument("workdir", help="Dataset directory (e.g., .../Px-N)")
    ap.add_argument(
        "--xna", "-x", required=True,
        help="One-letter XNA; BED at {workdir}/references/{XNA}.bed"
    )
    args = ap.parse_args()

    workdir = Path(args.workdir).resolve()
    xna = (args.xna or "").strip().upper()
    if len(xna) != 1 or not xna.isalpha():
        raise ValueError(f"[XNA] Must be a single letter A–Z. Got: {xna!r}")


    # Discover inputs
    can_root = workdir / "canonical_preprocess"
    mod_root = workdir / "modified_preprocess"
    can_pod5s = sorted((can_root / "pod5").glob("*.pod5"))
    mod_pod5s = sorted((mod_root / "pod5").glob("*.pod5"))
    if len(can_pod5s) != 1:
        details = "\n".join(f"  - {p}" for p in can_pod5s) or "  (none)"
        raise RuntimeError(f"[CAN_POD5] Expected exactly 1 .pod5, found {len(can_pod5s)} in {can_root/'pod5'}:\n{details}")
    if len(mod_pod5s) != 1:
        details = "\n".join(f"  - {p}" for p in mod_pod5s) or "  (none)"
        raise RuntimeError(f"[MOD_POD5] Expected exactly 1 .pod5, found {len(mod_pod5s)} in {mod_root/'pod5'}:\n{details}")

    can_pod5 = can_pod5s[0]
    mod_pod5 = mod_pod5s[0]
    can_bam  = can_root / "bam" / "aligned.sorted.bam"
    mod_bam  = mod_root / "bam" / "aligned.sorted.bam"
    xna_bed  = workdir / "references" / f"{xna}.bed"
    out_dir  = workdir / "signal_plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    for path, label in [(can_pod5, "CAN_POD5"), (mod_pod5, "MOD_POD5"), (can_bam, "CAN_BAM"), (mod_bam, "MOD_BAM"), (xna_bed, "XNA_BED")]:
        if not Path(path).exists():
            raise FileNotFoundError(f"[{label}] Expected file not found: {path}")


    # Window around XNA
    sites = load_xna_sites(str(xna_bed), CHROM)
    if not sites:
        raise ValueError(f"No XNA sites found for {CHROM} in {xna_bed}")
    xna_start, xna_end, _ = sites[0]
    start = xna_start - FLANK
    end   = xna_end + FLANK

    # Refinement per Remora notebooks
    sig_map_refiner = refine_signal_map.SigMapRefiner(
        kmer_model_filename=LEVELS_TXT,
        do_rough_rescale=True,
        scale_iters=0,
        do_fix_guage=True,
    )

    can_lines, kept_can, total_can = collect_lines(str(can_pod5), str(can_bam), CHROM, start, end, N_READS, sig_map_refiner)
    mod_lines, kept_mod, total_mod = collect_lines(str(mod_pod5), str(mod_bam), CHROM, start, end, N_READS, sig_map_refiner)
    n_can_pts = int(np.sum([ln[0].size for ln in can_lines])) if can_lines else 0
    n_mod_pts = int(np.sum([ln[0].size for ln in mod_lines])) if mod_lines else 0
    print(f"[DEBUG] Standard reads kept: {kept_can}/{total_can}  | total points: {n_can_pts}")
    print(f"[DEBUG] Modified reads kept: {kept_mod}/{total_mod} | total points: {n_mod_pts}")

    if not can_lines and not mod_lines:
        raise RuntimeError("No signal mapped in the requested region.")

    # Robust y-lims from data
    all_y = np.concatenate([ln[1] for ln in (can_lines + mod_lines)])
    lo = np.nanquantile(all_y, 0.01); hi = np.nanquantile(all_y, 0.99)
    pad = 0.06 * (hi - lo if hi > lo else 1.0)
    #ylo, yhi = lo - pad, hi + pad
    ylo, yhi = -2, 2

    # Figure
    fig, ax = plt.subplots(figsize=(5.9, 3.6))

    # Draw
    if RENDER_MODE == "line":
        for xs, ys in can_lines:
            if xs.size:
                ax.plot(xs, ys, color=COL_STD, alpha=ALPHA_STD, linewidth=0.75, rasterized=True)
        for xs, ys in mod_lines:
            if xs.size:
                ax.plot(xs, ys, color=COL_MOD, alpha=ALPHA_MOD, linewidth=0.75, rasterized=True)
    elif RENDER_MODE == "dots":
        for xs, ys in can_lines:
            if xs.size:
                ax.scatter(xs, ys, s=POINT_SIZE, alpha=ALPHA_STD, color=COL_STD, linewidths=0, rasterized=True)
        for xs, ys in mod_lines:
            if xs.size:
                ax.scatter(xs, ys, s=POINT_SIZE, alpha=ALPHA_MOD, color=COL_MOD, linewidths=0, rasterized=True)
    elif RENDER_MODE == "micro":
        for xs, ys in can_lines:
            if xs.size:
                micro_segments(ax, xs, ys, COL_STD, alpha=ALPHA_STD, seg_frac=0.65, lw=0.35)
        for xs, ys in mod_lines:
            if xs.size:
                micro_segments(ax, xs, ys, COL_MOD, alpha=ALPHA_MOD, seg_frac=0.65, lw=0.35)

    # XNA highlight (exact base)
    ax.axvspan(xna_start, xna_start + 1, color=COL_X_HIGHLIGHT, alpha=0.23, linewidth=0)
    if SHOW_XNA_CENTER_LINE:
        ax.axvline(xna_start + 0.5, color=COL_X_HIGHLIGHT, linestyle="--", linewidth=1.5)

    # Axes
    ax.set_xlim(start, end)
    ax.margins(x=0)        # removes default horizontal padding
    ax.set_ylim(ylo-0.4, yhi+0.4)   # robust range
    set_standard_signal_yticks(ax)

    set_base_center_ticks(ax, start, end, label_every=LABEL_EVERY)
    style_axes(ax, ylabel="Normalized signal")





    # Save
    stem = f"{workdir.name}_{xna}_{CHROM}_{start}-{end}_{RENDER_MODE}"
    svg = out_dir / f"{stem}.svg"
    pdf = out_dir / f"{stem}.pdf"
    plt.tight_layout(pad=0.9)
    plt.savefig(svg)
    plt.savefig(pdf)
    print(f"[INFO] Wrote figures: {svg} and {pdf}")

if __name__ == "__main__":
    main()

