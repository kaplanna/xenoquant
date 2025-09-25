#!/usr/bin/env python3
# xr_step_metric_pub.py
# Per-base step plot (median or trimmed mean) in reference coordinates, strand-aware.

from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pod5
import pysam
from remora import io, refine_signal_map

# ================================
# User-tunable params (paths resolved in main)
# ================================
# CLI example:
#   python xr_step_metric_pub.py /path/to/Px-N -x B --metric median
WORKDIR_DEFAULT = Path("/home/marchandlab/DataAnalysis/Kaplan/training/2509_Signal_Plots/250922_CNT_plots/Px-N")
XNA_DEFAULT     = "S"  # one-letter; BED at {WORKDIR}/references/{XNA}.bed

LEVELS_TXT = "/home/marchandlab/github/kaplanna/xemora/models/remora/9mer_10-4-1.tsv"
CHROM = "contig1"
N_READS = 50
FLANK = 10

# ---- Metric toggle ----
# Choose "median" or "trimmean"
METRIC  = "median"   # can be overridden via --metric
ST_TRIM = 1          # used only for trimmean
EN_TRIM = 1          # used only for trimmean

# Colors / styling
COL_STD = "#000000"   # Standard
COL_MOD = "#56B4E9"   # Modified
COL_XNA = "#9E9E9E"
SHOW_X_LABEL = False
SHOW_Y_LABEL = True
LABEL_EVERY  = None    # auto-thin x labels
LINE_W = 1.4
ALPHA = 0.95

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
    "savefig.transparent": True,
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
    ax.set_ylabel(ylabel if SHOW_Y_LABEL else "")

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

def per_read_metric(sm_region, ref_start, ref_end, strand, metric="median", st_trim=1, en_trim=1):
    """
    Compute per-base metric for one read over [ref_start, ref_end):
      - metric="median": per-base median of samples in the base slice
      - metric="trimmean": per-base trimmed mean (drop st_trim / en_trim samples if available)
    Returns array of length (ref_end - ref_start) aligned to genomic coordinates.
    """
    sig = sm_region.norm_signal
    s2s = sm_region.seq_to_sig_map
    nb  = len(sm_region.seq)
    n_sig = len(sig)
    L = ref_end - ref_start
    out = np.full(L, np.nan, dtype=float)

    def base_pos(b):
        return (ref_start + b) if strand == "+" else ((ref_end - 1) - b)

    for b in range(nb):
        s0 = s2s[b]
        s1 = s2s[b+1] if (b+1) < len(s2s) else None
        if s0 is None or s1 is None or s0 < 0 or s1 > n_sig or s1 <= s0:
            continue
        seg = sig[s0:s1]
        m = seg.size
        if m == 0:
            val = np.nan
        elif metric == "median":
            val = float(np.median(seg))
        elif metric == "trimmean":
            if m <= (st_trim + en_trim):
                val = float(np.mean(seg))
            else:
                val = float(np.mean(seg[st_trim: m - en_trim]))
        else:
            raise ValueError("metric must be 'median' or 'trimmean'")
        idx = base_pos(b) - ref_start
        if 0 <= idx < L:
            out[idx] = val
    return out

def collect_metric_arrays(pod5_path, bam_path, chrom, start, end, n_reads, refiner, metric="median", st_trim=1, en_trim=1):
    """
    Iterate reads overlapping [start, end), refine mapping, compute per-base metric per read.
    Returns stacked array (kept_reads, region_len).
    """
    dr = pod5.DatasetReader(pod5_path)
    bam_fh = pysam.AlignmentFile(bam_path, "rb")

    per_read = []
    kept = total = 0

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
        arr = per_read_metric(sm_region, start, end, strand, metric=metric, st_trim=st_trim, en_trim=en_trim)

        if not np.any(np.isfinite(arr)):
            continue

        per_read.append(arr)
        kept += 1
        if kept >= n_reads:
            break

    if kept == 0:
        return np.empty((0, end - start), dtype=float), kept, total
    return np.vstack(per_read), kept, total

def main():
    import argparse

    ap = argparse.ArgumentParser(description="Per-base step plot of normalized signal around an XNA site.")
    ap.add_argument("workdir", nargs="?", default=str(WORKDIR_DEFAULT), help="Dataset directory (e.g., .../Px-N)")
    ap.add_argument("--xna", "-x", default=XNA_DEFAULT, help="One-letter XNA; BED at {workdir}/references/{XNA}.bed")
    ap.add_argument("--metric", choices=("median", "trimmean"), default=METRIC, help="Per-base aggregation metric")
    ap.add_argument("--st-trim", type=int, default=ST_TRIM, help="Leading samples to trim (trimmean)")
    ap.add_argument("--en-trim", type=int, default=EN_TRIM, help="Trailing samples to trim (trimmean)")
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
    xna_bed  = workdir / "references" / f"{xna}.bed"  # <-- corrected folder
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
    L     = end - start

    # Refinement
    sig_map_refiner = refine_signal_map.SigMapRefiner(
        kmer_model_filename=LEVELS_TXT,
        do_rough_rescale=True,
        scale_iters=0,
        do_fix_guage=True,
    )

    can_arr, kept_can, total_can = collect_metric_arrays(
        str(can_pod5), str(can_bam), CHROM, start, end, N_READS, sig_map_refiner,
        metric=args.metric, st_trim=args.st_trim, en_trim=args.en_trim
    )
    mod_arr, kept_mod, total_mod = collect_metric_arrays(
        str(mod_pod5), str(mod_bam), CHROM, start, end, N_READS, sig_map_refiner,
        metric=args.metric, st_trim=args.st_trim, en_trim=args.en_trim
    )

    print(f"[DEBUG] Standard reads kept: {kept_can}/{total_can}, array shape: {can_arr.shape}")
    print(f"[DEBUG] Modified reads kept: {kept_mod}/{total_mod}, array shape: {mod_arr.shape}")

    if can_arr.size == 0 and mod_arr.size == 0:
        raise RuntimeError("No valid per-base metrics were computed in the region.")

    # Average across reads at each base (ignore NaNs)
    can_mean = np.nanmean(can_arr, axis=0) if can_arr.size else np.full(L, np.nan)
    mod_mean = np.nanmean(mod_arr, axis=0) if mod_arr.size else np.full(L, np.nan)

    # Robust y-limits
    pooled = np.concatenate([v[np.isfinite(v)] for v in (can_mean, mod_mean) if v.size])
    lo = np.quantile(pooled, 0.01) if pooled.size else -2
    hi = np.quantile(pooled, 0.99) if pooled.size else  2
    pad = 0.06 * (hi - lo if hi > lo else 1.0)
    ylo, yhi = lo - pad, hi + pad

    # Step edges
    def to_edges(vals):
        positions = np.arange(start, end)
        edges = np.concatenate([positions, [end]])
        v = np.concatenate([vals, [vals[-1] if vals.size else np.nan]])
        return edges, v

    edges_c, vals_c = to_edges(can_mean)
    edges_m, vals_m = to_edges(mod_mean)

    # Plot
    fig, ax = plt.subplots(figsize=(5.9, 3.6))
    ax.step(edges_c, vals_c, where="post", color=COL_STD, linewidth=LINE_W, alpha=ALPHA, label="Standard")
    ax.step(edges_m, vals_m, where="post", color=COL_MOD, linewidth=LINE_W, alpha=ALPHA, label="Modified")

    # XNA highlight
    ax.axvspan(xna_start, xna_start + 1, color=COL_XNA, alpha=0.23, linewidth=0)
    ax.axvline(xna_start + 0.5, color=COL_XNA, linestyle="--", linewidth=0.8)

    # Axes & ticks
    ax.set_xlim(start, end)
    ax.set_ylim(ylo, yhi)
    set_base_center_ticks(ax, start, end, label_every=LABEL_EVERY)
    style_axes(ax, xlabel=f"Position on {CHROM}", ylabel="Normalized signal")

    # Legend (top-right)
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=COL_STD, linewidth=LINE_W),
               Line2D([], [], color=COL_MOD, linewidth=LINE_W)]
    ax.legend(handles=handles, labels=["Standard", "Modified"],
              loc="upper right", fontsize=9.8, handlelength=2.0)

    # Save
    metric_tag = args.metric
    stem = f"{workdir.name}_{xna}_{CHROM}_{start}-{end}_step_{metric_tag}"
    svg = out_dir / f"{stem}.svg"
    pdf = out_dir / f"{stem}.pdf"
    plt.tight_layout(pad=0.9)
    plt.savefig(svg)
    plt.savefig(pdf)
#    plt.show()
    print(f"[INFO] Wrote figures: {svg} and {pdf}")

if __name__ == "__main__":
    main()

