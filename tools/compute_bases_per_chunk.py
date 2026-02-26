#!/usr/bin/env python3
# compute_bases_per_chunk.py
import logging
import pod5
import pysam
import numpy as np
from tqdm import tqdm

from remora import io, refine_signal_map

logging.getLogger("Remora").setLevel(logging.INFO)

# ==========================
# CONFIG — EDIT
# ==========================
POD5_PATH = "/home/marchandlab/DataAnalysis/Kaplan/basecall/PZn/251013_PG_Model_Testing/Can/PG_150_150_not-mixed/preprocess/pod5/50_pod5_test.pod5"
BAM_PATH  = "/home/marchandlab/DataAnalysis/Kaplan/basecall/PZn/251013_PG_Model_Testing/Can/PG_150_150_not-mixed/preprocess/bam/aligned.BAM"
BED_PATH  = "/home/marchandlab/DataAnalysis/Kaplan/basecall/PZn/251013_PG_Model_Testing/Can/PG_150_150_not-mixed/references/P.bed"
LEVELS    = "/home/marchandlab/github/kaplanna/xemora/models/remora/9mer_10-4-1.tsv"

K = 50            # number of signal samples upstream/downstream
MAX_READS = 500      # limit number of reads used

# ==========================
# HELPERS
# ==========================
def load_bed(path):
    sites = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            chrom, start, end, *rest = line.strip().split("\t")
            name = rest[0] if rest else f"{chrom}:{start}-{end}"
            sites.append((chrom, int(start), int(end), name))
    return sites

def get_pod5_record(reader: pod5.Reader, read_id: str):
    it = reader.reads(selection=[read_id], missing_ok=True, preload={"samples"})
    return next(it, None)

# ==========================
# MAIN
# ==========================
def main():
    sites = load_bed(BED_PATH)
    print(f"[STATUS] Loaded {len(sites)} XNA sites from {BED_PATH}")

    refiner = refine_signal_map.SigMapRefiner(
        kmer_model_filename=LEVELS,
        do_rough_rescale=True,
        scale_iters=0,
        do_fix_guage=True,
    )

    upstream_counts, downstream_counts = [], []
    seen_reads = set()
    n_sites_used = 0

    with pod5.Reader(POD5_PATH) as p5, pysam.AlignmentFile(BAM_PATH, "rb") as bam:
        with tqdm(total=MAX_READS, desc="Processing usable sites", unit="site") as pbar:
            for aln in bam.fetch(until_eof=True):
                if MAX_READS and len(seen_reads) >= MAX_READS:
                    break
                if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                    continue
                rid = aln.query_name
                if rid in seen_reads:
                    continue

                rec = get_pod5_record(p5, rid)
                if rec is None:
                    continue

                try:
                    io_read = io.Read.from_pod5_and_alignment(rec, aln)
                    io_read.set_refine_signal_mapping(refiner, ref_mapping=False)
                    io_read.set_refine_signal_mapping(refiner, ref_mapping=True)
                except Exception:
                    continue

                for chrom, start, end, name in sites:
                    if aln.reference_name != chrom:
                        continue
                    if not (aln.reference_start <= start < aln.reference_end):
                        continue

                    # Map BED site to alignment-relative base index
                    base_idx = start - aln.reference_start
                    try:
                        sm_reg = io_read.extract_ref_reg(io.RefRegion(
                            ctg=chrom,
                            strand=io_read.ref_reg.strand,
                            start=aln.reference_start,
                            end=aln.reference_end
                        ))
                    except Exception:
                        continue

                    s2s = sm_reg.seq_to_sig_map
                    if base_idx < 0 or base_idx + 1 >= len(s2s):
                        continue
                    s0, s1 = s2s[base_idx], s2s[base_idx + 1]
                    if s0 is None or s1 is None:
                        continue
                    center = (s0 + s1) // 2

                    # Define upstream and downstream signal windows
                    up_region = (center - K, center)
                    dn_region = (center, center + K)

                    bases_up, bases_dn = 0, 0
                    for b in range(len(sm_reg.seq)):
                        seg0, seg1 = s2s[b], s2s[b + 1] if b + 1 < len(s2s) else None
                        if seg0 is None or seg1 is None:
                            continue
                        if seg1 <= up_region[1] and seg0 >= up_region[0]:
                            bases_up += 1
                        if seg0 >= dn_region[0] and seg1 <= dn_region[1]:
                            bases_dn += 1

                    if bases_up > 0 or bases_dn > 0:
                        upstream_counts.append(bases_up)
                        downstream_counts.append(bases_dn)
                        n_sites_used += 1
                        pbar.update(1)

                seen_reads.add(rid)

    if not upstream_counts:
        print("[ERROR] No sites processed.")
        return

    up_arr, dn_arr = np.array(upstream_counts), np.array(downstream_counts)
    print(f"\n[RESULTS] Processed {n_sites_used} XNA sites from {len(seen_reads)} reads")
    print("Upstream bases covered:")
    print(f" median={np.median(up_arr):.2f}  IQR=[{np.percentile(up_arr,25):.2f}, {np.percentile(up_arr,75):.2f}]")
    print("Downstream bases covered:")
    print(f" median={np.median(dn_arr):.2f}  IQR=[{np.percentile(dn_arr,25):.2f}, {np.percentile(dn_arr,75):.2f}]")

if __name__ == "__main__":
    main()

