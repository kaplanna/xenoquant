#!/usr/bin/env python3
# xr_signal_metrics_discovery_wdorado.py
# ----------------------------------------------------------------------
# Genome-wide per-base signal metrics comparison (canonical vs modified)
# Automatically performs Dorado basecalling and Minimap2 alignment
# (with --emit-moves and --no-trim) if BAMs are missing or empty.
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
import subprocess, os, csv

# ================================================================
# === USER CONFIGURATION =========================================
# ================================================================

# Input POD5s
CAN_POD5 = Path("/home/marchandlab/DataAnalysis/Kaplan/training/2509_Signal_Plots/250928_GBC_Plots/BA-plots/canonical_preprocess/pod5/10_pod5.pod5")
MOD_POD5 = Path("/home/marchandlab/DataAnalysis/Kaplan/training/2509_Signal_Plots/250928_GBC_Plots/BA-plots/modified_preprocess/pod5/1_pod5.pod5")

# Reference FASTA
REF_FASTA = Path("/home/marchandlab/DataAnalysis/Kaplan/training/2509_Signal_Plots/250928_GBC_Plots/BA-plots/references/xGBC_90mer_clean_N.fa")

# Remora 9-mer k-mer level table
LEVELS_TXT = "/home/marchandlab/github/kaplanna/xemora/models/remora/9mer_10-4-1.tsv"

# Output directory
OUT_DIR = Path("/home/marchandlab/DataAnalysis/Kaplan/training/2509_Signal_Plots/"
               "250928_GBC_Plots/signal_discovery_testing")

# Dorado configuration (use your actual binary/model)
DORADO_PATH  = "~/dorado-0.8.0-linux-x64/bin/dorado"
DORADO_MODEL = "~/dorado-0.8.0-linux-x64/models/dna_r10.4.1_e8.2_400bps_hac@v5.0.0"
MIN_QSCORE   = 5
MAX_READS    = None

# Analysis settings
N_READS = 500
ST_TRIM, EN_TRIM = 1, 1
YLIMS = {
    "trimmean": (-2, 2),
    "diff": (-1, 1),
    "dwell": (0, 40),
    "trimsd": (0, 0.4),
}

# ================================================================
# === MATPLOTLIB STYLING ========================================
# ================================================================

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 11,
    "axes.labelsize": 11, "axes.titlesize": 11,
    "xtick.labelsize": 11, "ytick.labelsize": 11,
    "legend.fontsize": 11, "lines.linewidth": 2.0,
    "axes.linewidth": 2.0, "xtick.major.width": 2.0,
    "ytick.major.width": 2.0, "xtick.direction": "out",
    "ytick.direction": "out", "savefig.dpi": 300,
    "figure.dpi": 150, "pdf.fonttype": 42, "svg.fonttype": "none",
})
COL_STD, COL_MOD = "#4575b4", "#d73027"

# ================================================================
# === UTILITY FUNCTIONS ==========================================
# ================================================================

def run_cmd(cmd: str, soft_fail=False):
    """Run shell command and print output."""
    print(f"[CMD] {cmd}")
    rc = os.system(cmd)
    if rc != 0 and not soft_fail:
        raise RuntimeError(f"[ERROR] Command failed: {cmd}")
    return rc

def bam_is_empty(bam_path):
    """Return True if BAM missing or has no reads."""
    try:
        if not Path(bam_path).exists():
            return True
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for _ in bam.head(1):
                return False
        return True
    except Exception:
        return True

def dorado_basecall(
    dorado_path, dorado_model, pod5_file, bam_dir,
    min_qscore=None, max_reads=None, overwrite=False
):
    """Run Dorado basecaller → BAM with moves and no trimming."""
    dorado_path  = os.path.expanduser(dorado_path)
    dorado_model = os.path.expanduser(dorado_model)
    bam_dir = Path(bam_dir)
    bam_dir.mkdir(parents=True, exist_ok=True)
    out_bam = bam_dir / "basecalled.bam"

    if out_bam.exists() and not overwrite and not bam_is_empty(out_bam):
        print(f"[INFO] Skipping Dorado basecalling (BAM exists): {out_bam}")
        return out_bam

    args = f"--no-trim --emit-moves {pod5_file}"
    if max_reads:  args += f" -n {max_reads}"
    if min_qscore: args += f" --min-qscore {min_qscore}"

    cmd = f"{dorado_path} basecaller {dorado_model} {args} > {out_bam}"
    run_cmd(cmd)
    return out_bam

def minimap2_align_and_sort(bam_in, ref_fa, out_dir, regenerate=True):
    """Align Dorado BAM to reference, sort, and index."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    raw_bam, sorted_bam = out_dir / "aligned.unsorted.bam", out_dir / "aligned.sorted.bam"

    if sorted_bam.exists() and not regenerate and not bam_is_empty(sorted_bam):
        print(f"[INFO] Using existing sorted BAM: {sorted_bam}")
        return sorted_bam

    cmd_align = (
        f"samtools fastq -T '*' {bam_in} | "
        f"minimap2 -y -ax map-ont -k8 -w5 -s 80 --score-N 0 "
        f"--secondary no --sam-hit-only --MD {ref_fa} - | "
        f"samtools view -F0x800 -bho {raw_bam}"
    )
    run_cmd(cmd_align)
    run_cmd(f"samtools sort -o {sorted_bam} {raw_bam}")
    run_cmd(f"samtools index {sorted_bam}")
    if raw_bam.exists(): raw_bam.unlink()
    return sorted_bam

def ensure_index(bam_path):
    """Ensure BAM index exists."""
    if not (Path(str(bam_path) + ".bai").exists() or Path(str(bam_path) + ".csi").exists()):
        run_cmd(f"samtools index {bam_path}", soft_fail=True)

def parse_reference(fasta_path):
    record = next(SeqIO.parse(str(fasta_path), "fasta"))
    return record.id, len(record.seq)

def style_axes(ax, xlabel=None, ylabel=None):
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    if xlabel: ax.set_xlabel(xlabel)
    if ylabel: ax.set_ylabel(ylabel)

# ================================================================
# === SIGNAL PROCESSING HELPERS =================================
# ================================================================

def refine_and_extract(dr, aln, refiner, start, end):
    try:
        pod5_read = dr.get_read(aln.query_name)
        io_read = io.Read.from_pod5_and_alignment(pod5_read, aln)
        io_read.set_refine_signal_mapping(refiner, ref_mapping=False)
        io_read.set_refine_signal_mapping(refiner, ref_mapping=True)
        sm_region = io_read.extract_ref_reg(
            io_read.ref_reg.adjust(start_adjust=start - io_read.ref_reg.start,
                                   end_adjust=end - io_read.ref_reg.end))
        return sm_region, io_read.ref_reg.strand
    except Exception:
        return None, None

def per_base_metrics(sm_region, strand, region_len, st_trim=1, en_trim=1):
    sig, s2s = sm_region.norm_signal, sm_region.seq_to_sig_map
    nb = len(sm_region.seq)
    dwell, trimmean, trimsd = [], [], []
    for b in range(nb):
        s0, s1 = s2s[b], s2s[b + 1]
        if s0 is None or s1 is None or s1 <= s0: continue
        seg = sig[s0:s1]
        if seg.size == 0: continue
        seg_t = seg[st_trim: seg.size - en_trim] if seg.size > (st_trim + en_trim) else seg
        dwell.append(seg.size)
        trimmean.append(np.mean(seg_t))
        trimsd.append(np.std(seg_t, ddof=1) if seg_t.size > 1 else 0.0)
    arrs = [np.array(a) for a in (dwell, trimmean, trimsd)]
    if strand == "-": arrs = [a[::-1] for a in arrs]
    return {"dwell": arrs[0], "trimmean": arrs[1], "trimsd": arrs[2]}

def collect_metrics(pod5_path, bam_path, chrom, ref_len, refiner, n_reads):
    dr = pod5.DatasetReader(pod5_path)
    bam = pysam.AlignmentFile(bam_path, "rb")
    dwell, trimmean, trimsd, kept, total = [], [], [], 0, 0

    for aln in bam.fetch(chrom):
        total += 1
        sm_reg, strand = refine_and_extract(dr, aln, refiner, 0, ref_len)
        if sm_reg is None: continue
        met = per_base_metrics(sm_reg, strand, ref_len)
        for k in met:
            arr = met[k].astype(float)
            met[k] = np.pad(arr, (0, max(0, ref_len - len(arr))),
                            mode='constant', constant_values=np.nan)[:ref_len]
        dwell.append(met["dwell"])
        trimmean.append(met["trimmean"])
        trimsd.append(met["trimsd"])
        kept += 1
        if kept >= n_reads: break

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

    can_align_dir, mod_align_dir = OUT_DIR / "canonical_align", OUT_DIR / "modified_align"
    can_sorted_bam, mod_sorted_bam = can_align_dir / "aligned.sorted.bam", mod_align_dir / "aligned.sorted.bam"

    # --- Basecalling + alignment with auto-recovery ---
    need_can, need_mod = bam_is_empty(can_sorted_bam), bam_is_empty(mod_sorted_bam)
    if not need_can and not need_mod:
        print("[INFO] Both BAMs present and non-empty — skipping Dorado + Minimap2.")
    else:
        print("[INFO] Running Dorado basecalling + Minimap2 alignment as needed.")
        if need_can:
            print("[INFO] Canonical BAM missing/empty — regenerating.")
            can_bc_bam = dorado_basecall(DORADO_PATH, DORADO_MODEL, CAN_POD5,
                                         OUT_DIR / "canonical_basecall",
                                         min_qscore=MIN_QSCORE, max_reads=MAX_READS,
                                         overwrite=True)
            can_sorted_bam = minimap2_align_and_sort(can_bc_bam, REF_FASTA, can_align_dir)
        if need_mod:
            print("[INFO] Modified BAM missing/empty — regenerating.")
            mod_bc_bam = dorado_basecall(DORADO_PATH, DORADO_MODEL, MOD_POD5,
                                         OUT_DIR / "modified_basecall",
                                         min_qscore=MIN_QSCORE, max_reads=MAX_READS,
                                         overwrite=True)
            mod_sorted_bam = minimap2_align_and_sort(mod_bc_bam, REF_FASTA, mod_align_dir)

    ensure_index(can_sorted_bam)
    ensure_index(mod_sorted_bam)

    # --- Remora signal refinement + metric extraction ---
    refiner = refine_signal_map.SigMapRefiner(
        kmer_model_filename=LEVELS_TXT,
        do_rough_rescale=True, scale_iters=0, do_fix_guage=True,
    )

    can = collect_metrics(CAN_POD5, can_sorted_bam, chrom, ref_len, refiner, N_READS)
    mod = collect_metrics(MOD_POD5, mod_sorted_bam, chrom, ref_len, refiner, N_READS)

    mean_can = {k: np.nanmean(can[k], axis=0) for k in can}
    mean_mod = {k: np.nanmean(mod[k], axis=0) for k in mod}
    diff_trim = mean_mod["trimmean"] - mean_can["trimmean"]

    # --- Welch’s t-test ---
    pvals = np.full(ref_len, np.nan)
    for i in range(ref_len):
        cvals, mvals = can["trimmean"][:, i], mod["trimmean"][:, i]
        if np.sum(np.isfinite(cvals)) >= 2 and np.sum(np.isfinite(mvals)) >= 2:
            pvals[i] = ttest_ind(cvals[np.isfinite(cvals)],
                                 mvals[np.isfinite(mvals)],
                                 equal_var=False).pvalue
    neglog10p = -np.log10(pvals)

    # --- CSV output ---
    csv_out = OUT_DIR / f"{chrom}_signal_metrics.csv"
    with open(csv_out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["position", "mean_can", "mean_mod", "diff_trim",
                    "dwell_can", "dwell_mod", "sd_can", "sd_mod", "p_value"])
        for i in range(ref_len):
            w.writerow([
                i + 1,
                mean_can["trimmean"][i], mean_mod["trimmean"][i], diff_trim[i],
                mean_can["dwell"][i], mean_mod["dwell"][i],
                mean_can["trimsd"][i], mean_mod["trimsd"][i],
                pvals[i],
            ])
    print(f"[INFO] Wrote per-base metrics to {csv_out}")

    # --- Plot ---
    pos = np.arange(ref_len)
    fig, axes = plt.subplots(5, 1, figsize=(6, 10), sharex=True,
                             gridspec_kw={"hspace": 0.45})
    axA, axB, axC, axD, axE = axes


    from matplotlib.lines import Line2D
    legend_elems = [
        Line2D([0], [0], color=COL_STD, lw=2, label="Canonical"),
        Line2D([0], [0], color=COL_MOD, lw=2, label="Modified"),
    ]


    def step(ax, y, color):
        vals, edges = np.concatenate([y, [y[-1]]]), np.concatenate([pos, [ref_len]])
        ax.step(edges, vals, where="post", lw=2.0, color=color, alpha=0.95)

    step(axA, mean_can["trimmean"], COL_STD)
    step(axA, mean_mod["trimmean"], COL_MOD)
    axA.set_ylabel("Trimmed mean"); axA.set_ylim(*YLIMS["trimmean"])

    axB.axhline(0, color="#999999", ls="--", lw=1.5)
    step(axB, diff_trim, "#555555"); axB.set_ylabel("Δ Trimmed mean"); axB.set_ylim(*YLIMS["diff"])

    axC.scatter(pos+0.5, mean_can["dwell"], color=COL_STD, s=16)
    axC.scatter(pos+0.5, mean_mod["dwell"], color=COL_MOD, s=16)
    axC.set_ylabel("Dwell (samples)"); axC.set_ylim(*YLIMS["dwell"])

    axD.scatter(pos+0.5, mean_can["trimsd"], color=COL_STD, s=16)
    axD.scatter(pos+0.5, mean_mod["trimsd"], color=COL_MOD, s=16)
    axD.set_ylabel("Trimmed SD"); axD.set_ylim(*YLIMS["trimsd"])

    axE.scatter(pos+0.5, neglog10p, color="#444444", s=20)
    axE.axhline(-np.log10(0.05), color="#A9A9A9", ls="--", lw=1.5)
    axE.set_ylabel("−log10 p"); axE.set_xlabel("Reference position")

    for ax in axes:
        style_axes(ax); ax.margins(x=0); ax.set_xlim(0, ref_len)
    # After all subplots are created (before saving)
    fig.legend(
        handles=legend_elems,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.995),
        ncol=2,
        frameon=False,
    )


    pdf, svg = OUT_DIR / f"{chrom}_signal_metrics.pdf", OUT_DIR / f"{chrom}_signal_metrics.svg"
    plt.savefig(pdf); plt.savefig(svg)
    plt.close(fig)
    print(f"[INFO] Wrote plots:\n - {pdf}\n - {svg}")

# ================================================================
if __name__ == "__main__":
    main()

