#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import pod5
import pysam
from remora import io, refine_signal_map
import matplotlib as mpl

# ================================
# Config — EDIT THESE PATHS
# ================================
CAN_POD5 = "/home/marchandlab/DataAnalysis/Kaplan/training/DsPx/250922_CNT_plots/Ds-N/canonical_preprocess/pod5/pod5_train.pod5"
CAN_BAM  = "/home/marchandlab/DataAnalysis/Kaplan/training/DsPx/250922_CNT_plots/Ds-N/canonical_preprocess/bam/aligned.sorted.bam"
MOD_POD5 = "/home/marchandlab/DataAnalysis/Kaplan/training/DsPx/250922_CNT_plots/Ds-N/modified_preprocess/pod5/50_pod5_train.pod5"
MOD_BAM  = "/home/marchandlab/DataAnalysis/Kaplan/training/DsPx/250922_CNT_plots/Ds-N/modified_preprocess/bam/aligned.sorted.bam"
XNA_BED  = "/home/marchandlab/DataAnalysis/Kaplan/training/DsPx/250922_CNT_plots/Ds-N/references/B.bed"

# Still needed by Remora refiner, but NOT plotted anymore
LEVELS_TXT = "/home/marchandlab/github/kaplanna/xemora/models/remora/9mer_10-4-1.tsv"

CHROM = "contig1"
N_READS = 150
FLANK = 7
OUT_PREFIX = "/home/marchandlab/DataAnalysis/Kaplan/training/DsPx/250922_CNT_plots/Ds-N/ref_region_plots"
# ================================



def load_xna_sites(bed_file, chrom):
    sites = []
    with open(bed_file) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            fields = line.strip().split()
            if fields[0] == chrom:
                start = int(fields[1]); end = int(fields[2])
                name = fields[3] if len(fields) > 3 else "XNA"
                sites.append((start, end, name))
    return sites

import matplotlib as mpl

# Nature-ish, readable defaults
mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 10,           # bigger text
    "axes.labelsize": 11,
    "axes.titlesize": 11,
    "axes.linewidth": 0.8,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.size": 4,
    "ytick.major.size": 4,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "legend.frameon": False,
    "savefig.dpi": 300,
    "savefig.transparent": False,
    "svg.fonttype": "none",     # keep text as text
})


mpl.rcParams.update({
    # Keep text as text (not curves) in SVG
    "svg.fonttype": "none",    # ✅ editable SVG text
    # Keep TrueType in PDF so text stays editable
    "pdf.fonttype": 42,        # 42=TrueType, 3=Type3 (uneditable)
    "ps.fonttype": 42,         # for EPS, if you ever use it
    # Optional: avoid TeX text rendering; it often outlines glyphs
    "text.usetex": False,
})


# Okabe–Ito palette
COL_CAN = "#0072B2"   # Standard (blue)
COL_MOD = "#D55E00"   # Modified (vermillion)
COL_XNA = "#9E9E9E"   # XNA highlight gray

def style_axes(ax, *, xlabel=None, ylabel=None):
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_linewidth(0.8)
    ax.tick_params(axis="both", which="major", length=4, width=0.8, pad=3, labelsize=10)
    if xlabel is not None: ax.set_xlabel(xlabel, labelpad=3)
    if ylabel is not None: ax.set_ylabel(ylabel, labelpad=3)

def to_step_edges(x_positions, y_values):
    """
    Convert per-base values at integer positions x (len N) into step edges
    for where='post': edges len N+1 (last+1), values len N+1 (repeat last).
    """
    edges = np.concatenate([x_positions, [x_positions[-1] + 1]])
    y_step = np.concatenate([y_values, [y_values[-1]]])
    return edges, y_step
def set_base_end_ticks(ax, start, end, *, label_every=1):
    """
    Put ticks at base *edges* (integers): start, start+1, ... , end.
    Useful for non-step (spaghetti) plots.
    """
    edges = np.arange(start, end + 1)  # include right edge
    labels = [str(v) if (i % label_every) == 0 else "" for i, v in enumerate(edges)]
    ax.set_xticks(edges)
    ax.set_xticklabels(labels)

def set_base_center_ticks(ax, start, end, *, label_every=1):
    """
    Put a tick at every base center (i+0.5) from [start, end) and label with the
    reference coordinate i. If you need to thin labels, set label_every>1.
    """
    centers = np.arange(start, end) + 0.5
    labels = [str(p) if (k % label_every) == 0 else "" for k, p in enumerate(range(start, end))]
    ax.set_xticks(centers)
    ax.set_xticklabels(labels)

def per_position_valid_counts(traces):
    """Return x positions (ints) and count of finite values per base."""
    if not traces:
        return None, None
    X = traces[0][0]                              # integer ref positions [start, end)
    Y = np.stack([y for _, y, _ in traces], axis=0)
    counts = np.sum(np.isfinite(Y), axis=0)       # how many reads contributed per base
    return X, counts


        
        
        
def spaghetti_ref_region_dual(
    can_pod5, can_bam,
    mod_pod5, mod_bam,
    levels_txt,
    chrom, xna_start, xna_end,
    n_reads=200, flank=20,
    out_prefix="ref_region_plots"
):
    """Spaghetti plot + step plot (canonical vs modified means), no k-mer track."""

    # Readers
    can_pod5_dr = pod5.DatasetReader(can_pod5)
    mod_pod5_dr = pod5.DatasetReader(mod_pod5)
    can_bam_fh = pysam.AlignmentFile(can_bam, "rb")
    mod_bam_fh = pysam.AlignmentFile(mod_bam, "rb")

    # Refiner (uses provided model for mapping only)
    sig_map_refiner = refine_signal_map.SigMapRefiner(
        kmer_model_filename=levels_txt,
        do_rough_rescale=True,
        scale_iters=0,
        do_fix_guage=True,
    )

    # Define window around XNA
    start = xna_start - flank
    end = xna_end + flank
    expected_len = end - start

    # Tolerant collection: keep mostly-valid reads, mask bad bases with NaN
    def collect_traces(pod5_dr, bam_fh, max_reads, start, end,
                       max_invalid_frac=0., min_valid_steps=1):
        traces = []

        n_total = n_keyerr = n_region_short = n_bad_maplen = 0
        n_kept_strict = n_kept_tol = n_drop_invalid = 0

        for bam_aln in bam_fh.fetch(chrom, start, end):
            n_total += 1
            read_id = bam_aln.query_name
            try:
                pod5_read = pod5_dr.get_read(read_id)
            except KeyError:
                n_keyerr += 1
                continue

            io_read = io.Read.from_pod5_and_alignment(pod5_read, bam_aln)
            io_read.set_refine_signal_mapping(sig_map_refiner, ref_mapping=True)

            sm_region = io_read.extract_ref_reg(
                io_read.ref_reg.adjust(
                    start_adjust=start - io_read.ref_reg.start,
                    end_adjust=end   - io_read.ref_reg.end
                )
            )

            ref_seq = sm_region.seq
            if len(ref_seq) != expected_len:
                n_region_short += 1
                continue

            sig = sm_region.norm_signal
            sig_map = sm_region.seq_to_sig_map  # boundaries; expect len = len(seq)+1
            if len(sig_map) != (len(ref_seq) + 1):
                n_bad_maplen += 1
                continue

            base_signal = np.full(len(ref_seq), np.nan, dtype=float)
            invalid = 0
            for b in range(len(ref_seq)):
                s0 = sig_map[b]; s1 = sig_map[b + 1]
                if (s0 is None) or (s1 is None) or (s0 < 0) or (s1 > len(sig)) or (s1 <= s0):
                    invalid += 1
                else:
                    base_signal[b] = float(np.mean(sig[s0:s1]))

            frac_invalid = invalid / len(ref_seq)
            if invalid == 0:
                n_kept_strict += 1
                traces.append((np.arange(len(ref_seq)) + start, base_signal, ref_seq))
            elif (frac_invalid <= max_invalid_frac) and (1.0 - frac_invalid >= min_valid_steps):
                n_kept_tol += 1
                traces.append((np.arange(len(ref_seq)) + start, base_signal, ref_seq))
            else:
                n_drop_invalid += 1

            if len(traces) >= max_reads:
                break

        print(
            "[DEBUG] collect_traces (tolerant): "
            f"total={n_total}, kept={len(traces)} "
            f"(strict={n_kept_strict}, tolerant={n_kept_tol}), "
            f"keyerr={n_keyerr}, short={n_region_short}, bad_maplen={n_bad_maplen}, "
            f"dropped_for_invalid_frac={n_drop_invalid}"
        )
        return traces

    # Collect
    can_traces = collect_traces(can_pod5_dr, can_bam_fh, n_reads, start, end)
    mod_traces = collect_traces(mod_pod5_dr, mod_bam_fh, n_reads, start, end)
    print(f"[DEBUG] Collected {len(can_traces)} canonical traces")
    print(f"[DEBUG] Collected {len(mod_traces)} modified traces")

    if not can_traces and not mod_traces:
        raise RuntimeError("No valid traces found. Check alignments/mapping or relax tolerances.")

    # Dynamic y-limits from data (robust to outliers)
    def compute_ylim(traces_a, traces_b, q=0.005):
        ys = [tr[1] for tr in traces_a] + [tr[1] for tr in traces_b]
        if not ys:
            return (-2, 2)
        all_vals = np.concatenate(ys)
        all_vals = all_vals[np.isfinite(all_vals)]
        if all_vals.size == 0:
            return (-2, 2)
        lo = np.quantile(all_vals, q)
        hi = np.quantile(all_vals, 1 - q)
        pad = 0.05 * (hi - lo if hi > lo else 1.0)
        return (lo - pad, hi + pad)

    ylo, yhi = compute_ylim(can_traces, mod_traces)

    # ==========================
    # Plot 1: Spaghetti (pleasant 5.6×3.6 in)
    # ==========================
    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    for x, y, _ in can_traces:
        ax.plot(x, y, color=COL_CAN, alpha=0.14, linewidth=0.5)
    for x, y, _ in mod_traces:
        ax.plot(x, y, color=COL_MOD, alpha=0.14, linewidth=0.5)

    # Shade exactly the XNA base [start, start+1) and draw center line
    ax.axvspan(xna_start - 0.5, xna_start + 0.5, color=COL_XNA, alpha=0.30, linewidth=0)
    ax.axvline(xna_start, color=COL_XNA, linestyle="--", linewidth=0.8)

    # Legend (explicit handles/labels for MPL compatibility)
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=COL_CAN, linewidth=1.4),
               Line2D([], [], color=COL_MOD, linewidth=1.4)]
    labels = ["Standard", "Modified"]
    ax.legend(handles=handles, labels=labels, loc="upper left", fontsize=10, handlelength=2.2)

    # Limits: show the full base range as [start, end)
    ax.set_ylim(ylo, yhi)
    ax.set_xlim(start, end)

    # Ticks: one per base, centered at i
    set_base_end_ticks(ax, start, end, label_every=1)

    style_axes(ax, xlabel=f"Position on {chrom}", ylabel="Normalized signal")
    ax.set_title(f"{chrom}:{start}-{end}", pad=3)

    plt.tight_layout()
    plt.savefig(f"{out_prefix}_spaghetti.svg")
    plt.show()


    

    # ==========================
    # Plot 2: Step plot (means, 5.6×3.4 in)
    # ==========================
    def mean_trace(traces):
        if not traces: return None, None
        all_y = np.stack([y for _, y, _ in traces], axis=0)
        with np.errstate(all="ignore"):
            mean_y = np.nanmean(all_y, axis=0)
        x = traces[0][0]  # integer ref positions [start, end)
        return x, mean_y

    can_x, can_mean = mean_trace(can_traces)
    mod_x, mod_mean = mean_trace(mod_traces)

    fig, ax = plt.subplots(figsize=(5.6, 3.4))

    # Build edges so each base is the bin [i, i+1)
    if can_x is not None:
        can_edges, can_vals = to_step_edges(can_x, can_mean)
        ax.step(can_edges, can_vals, where="post", color=COL_CAN, linewidth=1.5, label="Standard")
    if mod_x is not None:
        mod_edges, mod_vals = to_step_edges(mod_x, mod_mean)
        ax.step(mod_edges, mod_vals, where="post", color=COL_MOD, linewidth=1.5, label="Modified")

    # XNA base
    ax.axvspan(xna_start, xna_start + 1, color=COL_XNA, alpha=0.25, linewidth=0)
    ax.axvline(xna_start + 0.5, color=COL_XNA, linestyle="--", linewidth=0.8)

    # Limits and ticks
    ax.set_ylim(ylo, yhi)
    ax.set_xlim(start, end)
    set_base_center_ticks(ax, start, end, label_every=1)

    ax.legend(loc="upper left", fontsize=10, handlelength=2.2)
    style_axes(ax, xlabel=f"Position on {chrom}", ylabel="Normalized signal")
    ax.set_title("Mean signal", pad=3)

    plt.tight_layout()
    plt.savefig(f"{out_prefix}_step.svg")
    plt.show()

    # ==========================
    # Plot 3: Number of reads per position
    # ==========================
    x_std, cnt_std = per_position_valid_counts(can_traces)
    x_mod, cnt_mod = per_position_valid_counts(mod_traces)

    fig, ax = plt.subplots(figsize=(5.6, 3.2))  # slightly wider than tall

    # Draw as step curves over base bins [i, i+1)
    if x_std is not None:
        e_std, v_std = to_step_edges(x_std, cnt_std)
        ax.step(e_std, v_std, where="post", color=COL_CAN, linewidth=1.5, label="Standard")
    if x_mod is not None:
        e_mod, v_mod = to_step_edges(x_mod, cnt_mod)
        ax.step(e_mod, v_mod, where="post", color=COL_MOD, linewidth=1.5, label="Modified")

    # Full base range and centered ticks at every base
    ax.set_xlim(start, end)
    set_base_center_ticks(ax, start, end, label_every=1)

    # Y limits: from 0 to max observed (pad a bit)
    max_cnt = 0
    if cnt_std is not None: max_cnt = max(max_cnt, int(np.nanmax(cnt_std)))
    if cnt_mod is not None: max_cnt = max(max_cnt, int(np.nanmax(cnt_mod)))
    ax.set_ylim(0, max(1, max_cnt) * 1.05)

    # XNA base highlight for orientation
    ax.axvspan(xna_start, xna_start + 1, color=COL_XNA, alpha=0.20, linewidth=0)
    ax.axvline(xna_start + 0.5, color=COL_XNA, linestyle="--", linewidth=0.8)

    ax.legend(loc="upper left", fontsize=10, handlelength=2.2)
    style_axes(ax, xlabel=f"Position on {chrom}", ylabel="# reads with signal")
    ax.set_title("Reads per position", pad=3)

    plt.tight_layout()
    plt.savefig(f"{out_prefix}_coverage.svg")
    plt.show()


    # ==========================
    # Plot N: Levels-style spaghetti (per-read flat steps)
    # ==========================
    fig, ax = plt.subplots(figsize=(5.6, 3.4))  # wider than tall, readable ratio

    # Draw each read as base-aligned steps [i, i+1)
    # Standard (canonical)
    for x, y, _ in can_traces:
        e, v = to_step_edges(x, y)
        ax.step(e, v, where="post", color=COL_CAN, alpha=0.10, linewidth=0.6)

    # Modified
    for x, y, _ in mod_traces:
        e, v = to_step_edges(x, y)
        ax.step(e, v, where="post", color=COL_MOD, alpha=0.10, linewidth=0.6)

    # XNA base highlight and center line (same as mean step plot)
    ax.axvspan(xna_start, xna_start + 1, color=COL_XNA, alpha=0.25, linewidth=0)
    ax.axvline(xna_start + 0.5, color=COL_XNA, linestyle="--", linewidth=0.8)

    # Limits and ticks match the step plot (ticks at centers)
    ax.set_ylim(ylo, yhi)
    ax.set_xlim(start, end)
    set_base_center_ticks(ax, start, end, label_every=1)

    # Minimal legend (optional)
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=COL_CAN, linewidth=1.2),
               Line2D([], [], color=COL_MOD, linewidth=1.2)]
    labels = ["Standard", "Modified"]
    ax.legend(handles=handles, labels=labels, loc="upper left", fontsize=10, handlelength=2.2)

    style_axes(ax, xlabel=f"Position on {chrom}", ylabel="Normalized signal")
    ax.set_title("Levels-style spaghetti (per read)", pad=3)

    plt.tight_layout()
    plt.savefig(f"{out_prefix}_levels_spaghetti.svg")
    plt.show()


if __name__ == "__main__":
    xna_sites = load_xna_sites(XNA_BED, CHROM)
    if not xna_sites:
        raise ValueError(f"No XNA sites found for {CHROM} in {XNA_BED}")
    xna_start, xna_end, xna_name = xna_sites[0]

    spaghetti_ref_region_dual(
        CAN_POD5, CAN_BAM,
        MOD_POD5, MOD_BAM,
        LEVELS_TXT,            # needed by Remora refiner only
        CHROM, xna_start, xna_end,
        n_reads=N_READS,
        flank=FLANK,
        out_prefix=OUT_PREFIX
    )

