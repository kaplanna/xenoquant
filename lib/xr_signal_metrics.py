#!/usr/bin/env python3
# xr_signal_metrics_fullcov.py
# Reference-anchored signal metrics with strict full-coverage per-read filtering,
# same refinement flow as spaghetti plot. No Remora batch metric calls.

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

# ================================
# User-tunable analysis params
# ================================
# You can override WORKDIR/XNA at the CLI:
#   python xr_signal_metrics_fullcov.py /path/to/Px-N -x S
WORKDIR_DEFAULT = Path("/home/marchandlab/DataAnalysis/Kaplan/training/2509_Signal_Plots/250922_CNT_plots/Px-N")
XNA_DEFAULT     = "S"  # one-letter base; BED must be {WORKDIR}/references/{XNA}.bed

LEVELS_TXT = "/home/marchandlab/github/kaplanna/xemora/models/remora/9mer_10-4-1.tsv"
CHROM = "contig1"
N_READS = 50
FLANK = 10

# Colors (Okabe–Ito)
COL_STD  = "#000000"  # Standard
COL_MOD  = "#56B4E9"  # Modified (sky blue)
COL_XNA  = "#9E9E9E"  # highlight fill
COL_DIFF = "#CC79A7"  # magenta (difference)

# Remora logging quieter
logging.getLogger("Remora").setLevel(logging.INFO)

# Matplotlib: minimalist, print-ready, editable text in SVG/PDF
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
    "svg.fonttype": "none",   # keep text as text for Illustrator
    "pdf.fonttype": 42,       # TrueType in PDF
    "ps.fonttype": 42,
    "text.usetex": False,
})

# ---------- small helpers (path-agnostic) ----------
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

def style_axes(ax, *, xlabel=None, ylabel=None):
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_linewidth(0.8)
    ax.tick_params(axis="both", which="major", length=3.8, width=0.8, pad=2.5, labelsize=10.5)
    ax.set_xlabel(xlabel if xlabel else "")
    ax.set_ylabel(ylabel if ylabel else "")

def set_base_center_ticks(ax, start, end):
    centers = np.arange(start, end) + 0.5
    width = end - start
    if width <= 35:      every = 1
    elif width <= 70:    every = 2
    elif width <= 150:   every = 5
    else:                every = 10
    labels = [str(p) if (k % every) == 0 else "" for k, p in enumerate(range(start, end))]
    ax.set_xticks(centers)
    ax.set_xticklabels(labels)

def refine_and_extract(pod5_dr, bam_aln, refiner, start, end):
    """Refine mapping and extract reference region. Return (sm_region, strand) or (None, None) on failure."""
    try:
        pod5_read = pod5_dr.get_read(bam_aln.query_name)
    except KeyError:
        return None, None
    try:
        io_read = io.Read.from_pod5_and_alignment(pod5_read, bam_aln)
        # parity with previous scripts & Remora notebook
        io_read.set_refine_signal_mapping(refiner, ref_mapping=False)
        io_read.set_refine_signal_mapping(refiner, ref_mapping=True)
        sm_region = io_read.extract_ref_reg(
            io_read.ref_reg.adjust(
                start_adjust=start - io_read.ref_reg.start,
                end_adjust=end   - io_read.ref_reg.end
            )
        )
        return sm_region, io_read.ref_reg.strand
    except Exception:
        return None, None

def full_coverage_metrics(sm_region, strand, expected_len, st_trim=1, en_trim=1):
    """
    Compute per-base metrics for a ReadRefReg with strict full coverage:
      - require len(seq)==expected_len and len(seq_to_sig_map)==len(seq)+1
      - require every interval valid (0 <= s0 < s1 <= len(sig))
    Flip arrays if strand=="-" to put them in increasing reference order.
    """
    sig = sm_region.norm_signal
    s2s = sm_region.seq_to_sig_map
    nb  = len(sm_region.seq)
    n_sig = len(sig)
    if nb != expected_len or len(s2s) != nb + 1:
        return None

    dwell   = np.empty(nb, dtype=float)
    mean    = np.empty(nb, dtype=float)
    trimmn  = np.empty(nb, dtype=float)
    trimsd  = np.empty(nb, dtype=float)

    for b in range(nb):
        s0 = s2s[b]; s1 = s2s[b + 1]
        if (s0 is None) or (s1 is None) or (s0 < 0) or (s1 > n_sig) or (s1 <= s0):
            return None
        seg = sig[s0:s1]
        dwell[b] = seg.size
        if seg.size == 0:
            return None
        mean[b] = float(np.mean(seg))
        if seg.size > (st_trim + en_trim):
            seg_t = seg[st_trim: seg.size - en_trim]
            trimmn[b] = float(np.mean(seg_t)) if seg_t.size > 0 else float(np.mean(seg))
            trimsd[b] = float(np.std(seg_t, ddof=1)) if seg_t.size > 1 else 0.0
        else:
            trimmn[b] = mean[b]
            trimsd[b] = float(np.std(seg, ddof=1)) if seg.size > 1 else 0.0

    if strand == "-":
        dwell  = dwell[::-1]
        mean   = mean[::-1]
        trimmn = trimmn[::-1]
        trimsd = trimsd[::-1]

    return {"dwell": dwell, "mean": mean, "trimmean": trimmn, "trimsd": trimsd}

def collect_sample_metrics(pod5_path, bam_path, chrom, start, end, n_reads, refiner):
    """Iterate BAM, refine & extract per-read; keep only reads with strict full coverage."""
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
        met = full_coverage_metrics(sm_reg, strand, expected_len=L, st_trim=1, en_trim=1)
        if met is None:
            continue
        for k in mats:
            mats[k].append(met[k])
        kept += 1
        if kept >= n_reads:
            break

    for k in mats:
        mats[k] = np.stack(mats[k], axis=0) if mats[k] else np.empty((0, L), dtype=float)

    print(f"[INFO] {Path(bam_path).stem}: kept={kept} (total overlapping={total}) with strict full coverage.")
    return mats

# ---------------- main (all path logic lives here) ----------------
def main():
    import argparse, sys

    ap = argparse.ArgumentParser(description="Signal metrics with strict full-coverage filtering.")
    ap.add_argument("workdir", nargs="?", default=str(WORKDIR_DEFAULT), help="Dataset directory (e.g., .../Px-N)")
    ap.add_argument("--xna", "-x", default=XNA_DEFAULT, help="One-letter XNA (BED at {workdir}/bed/{XNA}.bed)")
    args = ap.parse_args()

    # Resolve/validate inputs
    workdir = Path(args.workdir).resolve()
    xna = (args.xna or "").strip().upper()
    if len(xna) != 1 or not xna.isalpha():
        raise ValueError(f"[XNA] Must be a single letter A–Z. Got: {xna!r}")

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

    for path, label in [(can_pod5, "CAN_POD5"), (mod_pod5, "MOD_POD5"), (can_bam, "CAN_BAM"), (mod_bam, "MOD_BAM"), (xna_bed, "XNA_BED")]:
        if not Path(path).exists():
            raise FileNotFoundError(f"[{label}] Expected file not found: {path}")

    out_dir = workdir / "signal_plots"
    out_dir.mkdir(parents=True, exist_ok=True)

    print("[CONFIG] Using:")
    print(f"  WORKDIR : {workdir}")
    print(f"  XNA     : {xna}")
    print(f"  CAN_POD5: {can_pod5}")
    print(f"  CAN_BAM : {can_bam}")
    print(f"  MOD_POD5: {mod_pod5}")
    print(f"  MOD_BAM : {mod_bam}")
    print(f"  XNA_BED : {xna_bed}")
    print(f"  OUT_DIR : {out_dir}")

    # --- region from BED ---
    sites = load_xna_sites(str(xna_bed), CHROM)
    if not sites:
        raise ValueError(f"No XNA sites found for {CHROM} in {xna_bed}")
    xna_start, xna_end, xna_name = sites[0]
    start = xna_start - FLANK
    end   = xna_end + FLANK
    L     = end - start
    pos   = np.arange(start, end)

    # Naming stem for outputs
    out_stem = f"{workdir.name}_{xna}_{CHROM}_{start}-{end}"

    # --- Remora refiner ---
    sig_map_refiner = refine_signal_map.SigMapRefiner(
        kmer_model_filename=LEVELS_TXT,
        do_rough_rescale=True,
        scale_iters=0,
        do_fix_guage=True,
    )

    # --- metrics ---
    can = collect_sample_metrics(str(can_pod5), str(can_bam), CHROM, start, end, N_READS, sig_map_refiner)
    mod = collect_sample_metrics(str(mod_pod5), str(mod_bam), CHROM, start, end, N_READS, sig_map_refiner)

    # Coverage (should be flat at n_kept across bases)
    cov_can = np.sum(np.isfinite(can["trimmean"]), axis=0)
    cov_mod = np.sum(np.isfinite(mod["trimmean"]), axis=0)

    # Means across reads
    mean_can_trim = np.nanmean(can["trimmean"], axis=0)
    mean_mod_trim = np.nanmean(mod["trimmean"], axis=0)
    mean_can_dw   = np.nanmean(can["dwell"], axis=0)
    mean_mod_dw   = np.nanmean(mod["dwell"], axis=0)
    mean_can_sd   = np.nanmean(can["trimsd"], axis=0)
    mean_mod_sd   = np.nanmean(mod["trimsd"], axis=0)

    diff_trim = mean_mod_trim - mean_can_trim

    # Welch p-values per base on trimmed mean
    pvals = np.array([
        ttest_ind(can["trimmean"][:, i][np.isfinite(can["trimmean"][:, i])],
                  mod["trimmean"][:, i][np.isfinite(mod["trimmean"][:, i])],
                  equal_var=False).pvalue
        if (np.sum(np.isfinite(can["trimmean"][:, i])) >= 2 and
            np.sum(np.isfinite(mod["trimmean"][:, i])) >= 2)
        else np.nan
        for i in range(L)
    ], dtype=float)
    with np.errstate(divide="ignore"):
        neglog10p = -np.log10(pvals)

    # --- export TSV ---
    df = pl.DataFrame({
        "position": pos,
        "cov_standard": cov_can,
        "cov_modified": cov_mod,
        "trimmean_standard": mean_can_trim,
        "trimmean_modified": mean_mod_trim,
        "diff_trimmean": diff_trim,
        "dwell_standard": mean_can_dw,
        "dwell_modified": mean_mod_dw,
        "trimsd_standard": mean_can_sd,
        "trimsd_modified": mean_mod_sd,
        "neglog10p_trimmean": neglog10p,
    })
    tsv_path = out_dir / f"{out_stem}_remora_metrics.tsv"
    df.write_csv(str(tsv_path), separator="\t")
    print(f"[INFO] Wrote per-base metrics TSV: {tsv_path}")

    # --- plot ---
    fig_h = 7.4
    fig, axes = plt.subplots(
        nrows=5, ncols=1, figsize=(7.0, fig_h),
        gridspec_kw={"height_ratios": [1.2, 1.0, 1.0, 1.0, 1.0], "hspace": 0.28},
        sharex=True
    )
    axA, axB, axC, axD, axE = axes

    left_titles = [
        "A  Trimmed mean",
        "B  Difference (Mod − Std)",
        "C  Dwell",
        "D  Trimmed SD",
        "E  Coverage (reads/base)",
    ]
    for ax, title in zip(axes, left_titles):
        ax.text(-0.03, 1.02, title, transform=ax.transAxes,
                ha="left", va="bottom", fontsize=11, fontweight="bold")

    from matplotlib.lines import Line2D
    legend_handles = [
        Line2D([], [], color=COL_STD,  lw=1.6, label="Standard"),
        Line2D([], [], color=COL_MOD,  lw=1.6, label="Modified"),
        Line2D([], [], color=COL_DIFF, lw=1.6, label="Δ / -log10 p"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 1.02),
        fontsize=10
    )
    plt.subplots_adjust(top=0.92)

    edges = np.concatenate([pos, [end]])
    def step(ax, y, color, lw=1.2, alpha=0.95):
        vals = np.concatenate([y, [y[-1] if y.size else np.nan]])
        ax.step(edges, vals, where="post", linewidth=lw, alpha=alpha, color=color)

    # A) Trimmed mean (both)
    axA.axvspan(xna_start, xna_start + 1, color=COL_XNA, alpha=0.20, linewidth=0)
    step(axA, mean_can_trim, COL_STD)
    step(axA, mean_mod_trim, COL_MOD)
    axA.set_ylabel("Trimmed mean")

    # B) Difference
    axB.axvspan(xna_start, xna_start + 1, color=COL_XNA, alpha=0.20, linewidth=0)
    step(axB, diff_trim, COL_DIFF)
    axB.axhline(0.0, color="#666666", linewidth=0.8)
    axB.set_ylabel("Δ trim mean")

    # C) Dwell
    axC.axvspan(xna_start, xna_start + 1, color=COL_XNA, alpha=0.20, linewidth=0)
    step(axC, mean_can_dw, COL_STD)
    step(axC, mean_mod_dw, COL_MOD)
    axC.set_ylabel("Dwell")

    # D) Trimmed SD
    axD.axvspan(xna_start, xna_start + 1, color=COL_XNA, alpha=0.20, linewidth=0)
    step(axD, mean_can_sd, COL_STD)
    step(axD, mean_mod_sd, COL_MOD)
    axD.set_ylabel("Trim SD")

    # E) Coverage (lines). With strict full coverage, these should be flat.
    axE.axvspan(xna_start, xna_start + 1, color=COL_XNA, alpha=0.20, linewidth=0)
    step(axE, cov_can.astype(float), COL_STD, lw=1.2, alpha=0.9)
    step(axE, cov_mod.astype(float), COL_MOD, lw=1.2, alpha=0.9)
    axE.set_ylabel("Reads")

    for ax in axes:
        style_axes(ax)
    set_base_center_ticks(axE, start, end)
    axE.set_xlabel(f"Position on {CHROM}")

    plt.tight_layout(pad=0.9)
    svg = out_dir / f"{out_stem}_metrics.svg"
    pdf = out_dir / f"{out_stem}_metrics.pdf"
    plt.savefig(svg)
    plt.savefig(pdf)
#    plt.show()
    print(f"[INFO] Wrote figures: {svg} and {pdf}")

if __name__ == "__main__":
    main()

