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

# ================================
# User-tunable params (paths resolved in main)
# ================================
# You can override WORKDIR/XNA via CLI:
#   python xr_levels_spaghetti_lines_pub.py /path/to/Px-N -x B
WORKDIR_DEFAULT = Path("/home/marchandlab/DataAnalysis/Kaplan/training/2509_Signal_Plots/250922_CNT_plots/Px-N")
XNA_DEFAULT     = "S"  # one-letter; BED at {WORKDIR}/references/{XNA}.bed

LEVELS_TXT = "/home/marchandlab/github/kaplanna/xemora/models/remora/9mer_10-4-1.tsv"

CHROM = "contig1"
N_READS = 50
FLANK = 10

# Rendering mode: "line" (compressed-time polyline), "dots" (scatter), "micro" (tiny segments)
RENDER_MODE = "line"

# Colors / style
PALETTE_NAME = "black_sky"
PALETTES = {
    "blue_orange":    ("#0072B2", "#E69F00"),
    "green_magenta":  ("#009E73", "#CC79A7"),
    "sky_vermilion":  ("#56B4E9", "#D55E00"),
    "black_sky":      ("#000000", "#56b4e9"),
    "blue_vermilion": ("#0072B2", "#D55E00"),
}
COL_STD, COL_MOD = PALETTES[PALETTE_NAME]
COL_XNA = "#9E9E9E"    # highlight fill

# Plot tuning
POINT_SIZE  = 2.0
ALPHA_STD   = 0.15
ALPHA_MOD   = 0.15
BIN_JITTER  = 0.08     # 0..0.5 spreads samples within base bin
LABEL_EVERY = None     # auto if None
SHOW_X_LABEL = False
SHOW_Y_LABEL = False
SHOW_XNA_CENTER_LINE = False

# Minimal, print-ready MPL defaults
mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 10.5,
    "axes.labelsize": 11,
    "axes.titlesize": 11,
    "axes.linewidth": 0.8,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.size": 3.8,
    "ytick.major.size": 3.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "legend.frameon": False,
    "savefig.dpi": 300,
    "savefig.transparent": False,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "text.usetex": False,
})

def style_axes(ax, *, xlabel=None, ylabel=None):
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_linewidth(0.8)
    ax.tick_params(axis="both", which="major", length=3.8, width=0.8, pad=2.5, labelsize=10.5)
    ax.set_xlabel(xlabel if (xlabel and SHOW_X_LABEL) else "")
    ax.set_ylabel(ylabel if (ylabel and SHOW_Y_LABEL) else "")

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

    ap = argparse.ArgumentParser(description="Spaghetti lines of normalized signal around an XNA site.")
    ap.add_argument("workdir", nargs="?", default=str(WORKDIR_DEFAULT), help="Dataset directory (e.g., .../Px-N)")
    ap.add_argument("--xna", "-x", default=XNA_DEFAULT, help="One-letter XNA; BED at {workdir}/references/{XNA}.bed")
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

    print("[CONFIG] Using:")
    print(f"  WORKDIR : {workdir}")
    print(f"  XNA     : {xna}")
    print(f"  CAN_POD5: {can_pod5}")
    print(f"  CAN_BAM : {can_bam}")
    print(f"  MOD_POD5: {mod_pod5}")
    print(f"  MOD_BAM : {mod_bam}")
    print(f"  XNA_BED : {xna_bed}")
    print(f"  OUT_DIR : {out_dir}")

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
    ylo, yhi = lo - pad, hi + pad

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
    ax.axvspan(xna_start, xna_start + 1, color=COL_XNA, alpha=0.23, linewidth=0)
    if SHOW_XNA_CENTER_LINE:
        ax.axvline(xna_start + 0.5, color=COL_XNA, linestyle="--", linewidth=0.8)

    # Axes
    ax.set_xlim(start, end)
    ax.set_ylim(ylo, yhi)   # robust range
    set_base_center_ticks(ax, start, end, label_every=LABEL_EVERY)
    style_axes(ax, xlabel=f"Position on {CHROM}", ylabel="Normalized signal")

    # Legend
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=COL_STD, linewidth=1.2),
               Line2D([], [], color=COL_MOD, linewidth=1.2)]
    ax.legend(handles=handles, labels=["Standard", "Modified"],
              loc="upper right", fontsize=9.8, handlelength=2.0)

    # Save
    stem = f"{workdir.name}_{xna}_{CHROM}_{start}-{end}_{RENDER_MODE}"
    svg = out_dir / f"{stem}.svg"
    pdf = out_dir / f"{stem}.pdf"
    plt.tight_layout(pad=0.9)
    plt.savefig(svg)
    plt.savefig(pdf)
#    plt.show()
    print(f"[INFO] Wrote figures: {svg} and {pdf}")

if __name__ == "__main__":
    main()

