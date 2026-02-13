#!/usr/bin/env python3

import pysam
import os
import sys

# ============================
# USER SETTINGS
# ============================
INPUT_BAM = "/home/marchandlab/DataAnalysis/Kaplan/basecall/8letter/251013_P8_H4_Basecall/ZC-Basecall/preprocess/bam/aligned.BAM"


# Optional but STRONGLY recommended for IGV
INPUT_FASTA = "/home/marchandlab/DataAnalysis/Kaplan/basecall/8letter/251013_P8_H4_Basecall/PG-75_75_Basecall/references/xREF_P8_H4_PN.fasta"


OUTDIR = "/home/marchandlab/DataAnalysis/Kaplan/basecall/8letter/251013_P8_H4_Basecall/cleaned_bam_igv"
OUT_PREFIX = "cleaned_contigs"

# ============================
# Create output directory
# ============================
os.makedirs(OUTDIR, exist_ok=True)

# ============================
# Output file paths
# ============================
OUT_BAM = os.path.join(OUTDIR, f"{OUT_PREFIX}.bam")
OUT_FASTA = os.path.join(OUTDIR, f"{OUT_PREFIX}.fasta")
MAP_TSV = os.path.join(OUTDIR, f"{OUT_PREFIX}_contig_map.tsv")

# ============================
# Sanity checks
# ============================
if not os.path.exists(INPUT_BAM):
    sys.exit(f"[ERROR] BAM not found: {INPUT_BAM}")

if INPUT_FASTA is not None and not os.path.exists(INPUT_FASTA):
    sys.exit(f"[ERROR] FASTA not found: {INPUT_FASTA}")

# ============================
# Load BAM and header
# ============================
bam = pysam.AlignmentFile(INPUT_BAM, "rb")
header = bam.header.to_dict()

if "SQ" not in header:
    sys.exit("[ERROR] BAM header contains no @SQ entries")

# ============================
# Build contig mapping
# ============================
orig_contigs = [sq["SN"] for sq in header["SQ"]]

contig_map = {
    orig: f"contig{i+1}"
    for i, orig in enumerate(orig_contigs)
}

# ============================
# Save mapping
# ============================
with open(MAP_TSV, "w") as f:
    f.write("original_contig\tnew_contig\n")
    for orig, new in contig_map.items():
        f.write(f"{orig}\t{new}\n")

print(f"[INFO] Wrote contig map → {MAP_TSV}")

# ============================
# Rewrite BAM header
# ============================
for sq in header["SQ"]:
    old = sq["SN"]
    sq["SN"] = contig_map[old]

# ============================
# Write new BAM
# ============================
out_bam = pysam.AlignmentFile(OUT_BAM, "wb", header=header)

for read in bam:
    out_bam.write(read)

bam.close()
out_bam.close()

# ============================
# Sort BAM before indexing
# ============================
SORTED_BAM = os.path.join(OUTDIR, f"{OUT_PREFIX}.sorted.bam")

pysam.sort(
    "-o", SORTED_BAM,
    OUT_BAM
)

pysam.index(SORTED_BAM)

print(f"[INFO] Sorted BAM → {SORTED_BAM}")
print(f"[INFO] Indexed BAM")



# ============================
# Optional: rewrite FASTA
# ============================
if INPUT_FASTA is not None:
    print("[INFO] Rewriting FASTA to match contig names")

    fasta_in = pysam.FastaFile(INPUT_FASTA)

    with open(OUT_FASTA, "w") as out:
        for orig, new in contig_map.items():
            seq = fasta_in.fetch(orig)
            out.write(f">{new}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")

    fasta_in.close()

    pysam.faidx(OUT_FASTA)

    print(f"[INFO] Wrote cleaned FASTA → {OUT_FASTA}")
    print(f"[INFO] Indexed FASTA")

# ============================
# Done
# ============================
print("[DONE] Contig cleanup complete")
print("All outputs saved to:")
print(f"  {OUTDIR}")
print("\nLoad into IGV in this order:")
print("  1) FASTA (if provided)")
print("  2) BAM")


