#!/usr/bin/env python3
"""
Batch-merge two POD5 files across a range of ratios with replicates.

- Creates a directory for each ratio (e.g., fracB_0.30)
  and subdirectories rep1, rep2, rep3 with a merged .pod5 in each.
- Fractions are normalized automatically.
- Sampling is without replacement within each output; different outputs are independent.
- If one input lacks enough reads to hit the target, the other tops up (with warnings).

Requires:
    pip install pod5
"""

# -------------------------
# CONFIG — EDIT THESE
# -------------------------
POD5_A = "/home/marchandlab/DataAnalysis/Kaplan/raw/BS_xr/250207_GAC_PAGE_xr_train/20250207_1302_MN37138_AXE551_c0ddb701/pod5.pod5"
POD5_B = "/home/marchandlab/DataAnalysis/Kaplan/raw/BS_xr/250207_GBC_xr_train/20250207_1305_MN41475_AXE602_2b402855/pod5.pod5"

# Root directory to hold all ratio folders
OUTPUT_ROOT = "/home/marchandlab/DataAnalysis/Kaplan/basecall/Standard_Curve/mixed_pod5"

TOTAL_READS = 200_000        # target reads per output file
RATIO_STEP = 0.10            # step size for frac_B (e.g., 0.10 → 0.00, 0.10, ..., 1.00). Set to 1.0 for just 0 and 1.
REPLICATES = 3               # number of replicates per ratio

# Random seed: None → fully random each run (preferred for independent standard-curve points)
SEED = None

# -------------------------
# IMPLEMENTATION
# -------------------------
import os
import sys
import math
import random
from typing import List, Tuple, Set

import pod5 as p5  # pip install pod5


def get_read_ids(pod5_path: str) -> List[str]:
    """Return all read_ids (as strings) present in a POD5 file."""
    ids = []
    with p5.Reader(pod5_path) as reader:
        for rec in reader.reads():
            ids.append(str(rec.read_id))
    return ids


def sample_ids(ids: List[str], k: int, rng: random.Random) -> List[str]:
    """Sample k unique ids without replacement; if k >= len(ids), return all."""
    if k >= len(ids):
        return list(ids)
    return rng.sample(ids, k)


def compute_counts(
    n_total: int,
    frac_a: float,
    frac_b: float,
    count_a_avail: int,
    count_b_avail: int,
) -> Tuple[int, int, List[str]]:
    """
    Decide how many reads to draw from each file given desired fractions
    and availability. Fractions are normalized if they don't sum to 1.
    If one file is short, tries to top up from the other.
    """
    if frac_a <= 0 and frac_b <= 0:
        raise ValueError("At least one fraction must be > 0.")

    s = frac_a + frac_b
    frac_a_norm = frac_a / s
    frac_b_norm = frac_b / s

    # initial targets
    a_target = round(n_total * frac_a_norm)
    b_target = n_total - a_target

    warnings = []

    # clip by availability
    a_pick = min(a_target, count_a_avail)
    b_pick = min(b_target, count_b_avail)

    # top up if short
    shortfall = n_total - (a_pick + b_pick)
    if shortfall > 0:
        a_left = count_a_avail - a_pick
        from_a = min(shortfall, max(a_left, 0))
        a_pick += from_a
        shortfall -= from_a

        if shortfall > 0:
            b_left = count_b_avail - b_pick
            from_b = min(shortfall, max(b_left, 0))
            b_pick += from_b
            shortfall -= from_b

    if a_pick + b_pick < n_total:
        warnings.append(
            f"Requested total {n_total:,} but only {a_pick + b_pick:,} reads are available "
            f"({a_pick:,} from A, {b_pick:,} from B)."
        )
    if a_pick < a_target:
        warnings.append(f"File A limited: wanted {a_target:,}, taking {a_pick:,}.")
    if b_pick < b_target:
        warnings.append(f"File B limited: wanted {b_target:,}, taking {b_pick:,}.")

    return a_pick, b_pick, warnings


def write_reads(
    out_path: str,
    path_a: str,
    path_b: str,
    sel_a: Set[str],
    sel_b: Set[str],
    software_name: str = "merge_pod5_by_ratio",
) -> None:
    """
    Write selected reads from A and B into a single output POD5.

    Attempts fast selection using `Reader.reads(sel_ids)`. If that API
    isn’t supported in your installed pod5 version, falls back to scanning
    all reads and filtering in Python.
    """
    # ensure parent directory
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    with p5.Writer(out_path, software_name=software_name) as writer:
        # --- A ---
        if sel_a:
            try:
                with p5.Reader(path_a) as ra:
                    for rec in ra.reads(sel_a):
                        writer.add_read(rec.to_read())
            except TypeError:
                with p5.Reader(path_a) as ra:
                    for rec in ra.reads():
                        if str(rec.read_id) in sel_a:
                            writer.add_read(rec.to_read())

        # --- B ---
        if sel_b:
            try:
                with p5.Reader(path_b) as rb:
                    for rec in rb.reads(sel_b):
                        writer.add_read(rec.to_read())
            except TypeError:
                with p5.Reader(path_b) as rb:
                    for rec in rb.reads():
                        if str(rec.read_id) in sel_b:
                            writer.add_read(rec.to_read())


def frange_0_to_1(step: float) -> List[float]:
    """Generate [0.0, step, 2*step, ..., 1.0] with robust rounding."""
    if step <= 0 or step > 1:
        raise ValueError("RATIO_STEP must be in (0, 1].")
    n = int(round(1.0 / step))
    vals = [round(i * step, 8) for i in range(n + 1)]
    # ensure 1.0 exactly present
    vals[-1] = 1.0
    # ensure 0.0 exactly present
    vals[0] = 0.0
    return vals


def main() -> None:
    # Single RNG instance so subsequent draws are independent but not reproducible
    rng = random.Random() if SEED is None else random.Random(SEED)

    print("[INFO] Scanning input files once…")
    ids_a = get_read_ids(POD5_A)
    ids_b = get_read_ids(POD5_B)
    print(f"[INFO] File A: {len(ids_a):,} reads")
    print(f"[INFO] File B: {len(ids_b):,} reads")

    # Prepare ratios (iterate frac_B; frac_A = 1 - frac_B)
    ratios_b = frange_0_to_1(RATIO_STEP)

    for frac_b in ratios_b:
        frac_a = 1.0 - frac_b
        ratio_dir = os.path.join(OUTPUT_ROOT, f"fracB_{frac_b:0.2f}")
        os.makedirs(ratio_dir, exist_ok=True)

        # Decide intended counts for this ratio based on availability
        a_pick, b_pick, warns = compute_counts(
            n_total=TOTAL_READS,
            frac_a=frac_a,
            frac_b=frac_b,
            count_a_avail=len(ids_a),
            count_b_avail=len(ids_b),
        )

        # Log ratio header
        print("\n" + "=" * 80)
        print(f"[RATIO] B={frac_b:0.2f}  A={frac_a:0.2f}  → targets: A={a_pick:,}, B={b_pick:,}")
        for w in warns:
            print(f"[WARN] {w}", file=sys.stderr)

        for rep in range(1, REPLICATES + 1):
            rep_dir = os.path.join(ratio_dir, f"rep{rep}")
            os.makedirs(rep_dir, exist_ok=True)

            # Sample independently for each replicate
            sel_a = set(sample_ids(ids_a, a_pick, rng))
            sel_b = set(sample_ids(ids_b, b_pick, rng))

            # Avoid accidental overlap in read IDs (vanishingly unlikely)
            overlap = sel_a & sel_b
            if overlap:
                print(f"[WARN] {len(overlap)} duplicate read IDs across A and B; dropping from B.")
                sel_b -= overlap

            # Construct output filename
            pct_a = int(round(frac_a * 100))
            pct_b = int(round(frac_b * 100))
            total_tag = f"{TOTAL_READS//1000}k" if TOTAL_READS % 1000 == 0 else f"{TOTAL_READS}"
            out_name = f"mixed_A{pct_a}_B{pct_b}_rep{rep}_{total_tag}.pod5"
            out_path = os.path.join(rep_dir, out_name)

            actual_total = len(sel_a) + len(sel_b)
            print(f"[RUN ] B={frac_b:0.2f} rep{rep}: writing {actual_total:,} → {out_path}")
            write_reads(out_path, POD5_A, POD5_B, sel_a, sel_b)

    print("\n[DONE] All ratios and replicates completed.")
    print(f"[OUTPUT ROOT] {OUTPUT_ROOT}")


if __name__ == "__main__":
    main()

