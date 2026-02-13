#!/usr/bin/env python3
# xr_signal_extract_xna_base.py
# Extract per-read signal metrics at the exact XNA position (single base, no flanks)

from pathlib import Path
import numpy as np
import polars as pl
import pod5
import pysam
from remora import io, refine_signal_map
from xr_params import FLANK  # ensures BED compatibility if needed

# ============================
# User parameters
# ============================
WORKDIR_DEFAULT = Path("/home/marchandlab/DataAnalysis/Kaplan/training/2509_Signal_Plots/250922_CNT_plots/Px-N")
XNA_DEFAULT     = "S"
CHROM           = "contig1"
LEVELS_TXT      = "/home/marchandlab/github/kaplanna/xemora/models/remora/9mer_10-4-1.tsv"
N_READS         = 200
ST_TRIM         = 1
EN_TRIM         = 1

# ============================
# Helpers
# ============================
def load_xna_site(bed_file, chrom):
    with open(bed_file) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.strip().split()
            if fields[0] == chrom:
                start = int(fields[1]); end = int(fields[2])
                name = fields[3] if len(fields) > 3 else "XNA"
                return start, end, name
    raise ValueError(f"No XNA site found for {chrom} in {bed_file}")

def refine_and_extract(dr, aln, refiner, start, end):
    """Refine Remora mapping and extract region [start, end)."""
    try:
        pod5_read = dr.get_read(aln.query_name)
        io_read = io.Read.from_pod5_and_alignment(pod5_read, aln)
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

def get_single_base_metrics(sm_region, strand, base_idx, st_trim=1, en_trim=1):
    """Return dwell, mean, trimmed mean, and normalized signal for one base."""
    sig = sm_region.norm_signal
    s2s = sm_region.seq_to_sig_map
    nb = len(sm_region.seq)
    if base_idx < 0 or base_idx >= nb:
        return None

    s0, s1 = s2s[base_idx], s2s[base_idx + 1]
    if s0 is None or s1 is None or s1 <= s0:
        return None
    seg = sig[s0:s1]
    dwell = seg.size
    if dwell == 0:
        return None
    seg_t = seg[st_trim: seg.size - en_trim] if seg.size > (st_trim + en_trim) else seg
    mean = np.mean(seg)
    trimmean = np.mean(seg_t)
    # Return all signal values as comma-separated string for Excel histogram
    sig_str = ",".join(f"{v:.4f}" for v in seg.tolist())
    return dwell, mean, trimmean, sig_str

# ============================
# Main
# ============================
def main():
    import argparse
    ap = argparse.ArgumentParser(description="Extract signal metrics at XNA position per read.")
    ap.add_argument("workdir", nargs="?", default=str(WORKDIR_DEFAULT))
    ap.add_argument("--xna", "-x", default=XNA_DEFAULT)
    ap.add_argument("--nreads", "-n", type=int, default=N_READS)
    args = ap.parse_args()

    workdir = Path(args.workdir).resolve()
    xna = args.xna.upper()
    n_reads = args.nreads

    can_root = workdir / "canonical_preprocess"
    mod_root = workdir / "modified_preprocess"
    can_pod5 = next((can_root / "pod5").glob("*.pod5"))
    mod_pod5 = next((mod_root / "pod5").glob("*.pod5"))
    can_bam  = can_root / "bam" / "aligned.sorted.bam"
    mod_bam  = mod_root / "bam" / "aligned.sorted.bam"
    bed_file = workdir / "references" / f"{xna}.bed"

    out_dir = workdir / "signal_extract"
    out_dir.mkdir(parents=True, exist_ok=True)

    xna_start, xna_end, xname = load_xna_site(bed_file, CHROM)
    # For one-base site, use start as base index window
    start, end = xna_start, xna_end

    refiner = refine_signal_map.SigMapRefiner(
        kmer_model_filename=LEVELS_TXT,
        do_rough_rescale=True,
        scale_iters=0,
        do_fix_guage=True,
    )

    # ------------------------
    # Iterate both standard and modified
    # ------------------------
    records = []
    for label, pod5_path, bam_path in [("standard", can_pod5, can_bam), ("modified", mod_pod5, mod_bam)]:
        dr = pod5.DatasetReader(pod5_path)
        bam = pysam.AlignmentFile(bam_path, "rb")
        kept = 0
        for aln in bam.fetch(CHROM, start, end):
            if kept >= n_reads:
                break
            sm_reg, strand = refine_and_extract(dr, aln, refiner, start, end)
            if sm_reg is None:
                continue
            nb = len(sm_reg.seq)
            # Identify which base corresponds to the XNA site
            # (XNA position = first base in sm_region.seq if trimmed correctly)
            base_idx = 0 if strand == "+" else nb - 1
            metrics = get_single_base_metrics(sm_reg, strand, base_idx)
            if metrics is None:
                continue
            dwell, mean, trimmean, sig_str = metrics
            records.append({
                "read_id": aln.query_name,
                "label": label,
                "strand": strand,
                "dwell": dwell,
                "mean_signal": mean,
                "trimmed_mean": trimmean,
                "normalized_signal_values": sig_str,
            })
            kept += 1
        print(f"[INFO] {label}: extracted {kept} reads at XNA base")

    df = pl.DataFrame(records)
    out_path = out_dir / f"{workdir.name}_{xna}_{CHROM}_XNA_base_metrics.tsv"
    df.write_csv(out_path, separator="\t")
    print(f"[INFO] Wrote per-read base metrics: {out_path}")

if __name__ == "__main__":
    main()

