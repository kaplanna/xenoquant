#!/usr/bin/env python3
########################################################################
"""
xr_train_methods_bca.py  —  Xemora single-context training (Basecall-Anchored)

Hybrid workflow:
  • Uses BED file to define XNA site(s)
  • Extracts motif context from reference FASTA around each site
  • Performs basecall-anchored Remora chunk generation for both
    modified (X) and canonical (A/T/G/C) versions of that motif

All other stages (FASTA conversion, Dorado basecalling, alignment,
training) remain identical to the prior reference-anchored pipeline.
"""
########################################################################

import os
import sys
import re
import shlex
import subprocess
from pathlib import Path
from Bio import SeqIO
from typing import Optional
from xr_tools import *
from xr_params import *



print("Xemora [STATUS] - Initializing basecall-anchored training...")

# ---------------------------------------------------------------------
# CLI args
# ---------------------------------------------------------------------
if len(sys.argv) < 6:
    print("Usage: python xr_train_methods.py <workdir> <xna_raw_dir> <xna_ref_fasta> <dna_raw_dir> <dna_ref_fasta>")
    sys.exit(1)

working_dir   = os.path.expanduser(sys.argv[1])
xna_raw_dir   = os.path.expanduser(sys.argv[2])
xna_ref_fasta = os.path.expanduser(sys.argv[3])
dna_raw_dir   = os.path.expanduser(sys.argv[4])
dna_ref_fasta = os.path.expanduser(sys.argv[5])

# ---------------------------------------------------------------------
# Paths & directories
# ---------------------------------------------------------------------
working_dir = check_make_dir(working_dir)
ref_dir     = check_make_dir(os.path.join(working_dir, "references"))
chunk_dir   = check_make_dir(os.path.join(working_dir, "chunks"))
model_dir   = check_make_dir(os.path.join(working_dir, "model"))

mod_dir     = check_make_dir(os.path.join(working_dir, "modified_preprocess"))
mod_pod_dir = check_make_dir(os.path.join(mod_dir, "pod5"))
mod_bam_dir = check_make_dir(os.path.join(mod_dir, "bam"))

can_dir     = check_make_dir(os.path.join(working_dir, "canonical_preprocess"))
can_pod_dir = check_make_dir(os.path.join(can_dir, "pod5"))
can_bam_dir = check_make_dir(os.path.join(can_dir, "bam"))

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------
# Basic utilities
# ---------------------------------------------------------------------
def resolve_param_path(p: str) -> str:
    return os.path.abspath(os.path.join(BASE_DIR, os.path.expanduser(p)))

def validate_read_directory(raw_dir):
    if not os.path.isdir(raw_dir):
        print(f"[ERROR] Directory not found: {raw_dir}")
        sys.exit(1)
    exts = {f.split(".")[-1] for f in os.listdir(raw_dir) if "." in f}
    if len(exts) != 1:
        print(f"[WARN] Mixed filetypes in {raw_dir}: {sorted(exts)}")
        return ""
    return list(exts)[0]

def cod5_to_fast5(fast5_input, pod_dir, overwrite_pod):
    out = os.path.join(pod_dir, os.path.basename(fast5_input) + ".pod5")
    if overwrite_pod or not os.path.exists(out):
        os.system(f"pod5 convert fast5 --force-overwrite {fast5_input}/*.fast5 -o {out}")
    return out

def pod5_merge(pod5_input, pod_dir, overwrite_pod):
    out = os.path.join(pod_dir, os.path.basename(pod5_input) + ".pod5")
    if overwrite_pod or not os.path.exists(out):
        os.system(f"pod5 merge --force-overwrite {pod5_input}/*.pod5 -o {out}")
    return out

def dorado_basecall(dorado_path, dorado_model, min_qscore, pod5_file, bam_dir, basecall_pod, max_reads, filter_ids):
    out_bam = os.path.join(bam_dir, "bc.bam")
    if basecall_pod or not os.path.exists(out_bam):
        print("Xemora [STATUS] - Performing Dorado basecalling")
        args = f"--no-trim --emit-moves {pod5_file}"
        if filter_ids: args += f" -l {filter_ids}"
        if max_reads: args += f" -n {max_reads}"
        if min_qscore: args += f" --min-qscore {min_qscore}"
        os.system(f"{dorado_path} basecaller {dorado_model} {args} > {out_bam}")
    else:
        print("Xemora [STATUS] - Skipping Dorado basecall")
    return out_bam

def minimap2_aligner(input_bam, xfasta_path, bam_directory):
    sorted_bam = os.path.join(bam_directory, "aligned.sorted.bam")
    if regenerate_bam:
        cmd = (
            f"samtools fastq -T '*' {input_bam} | "
            f"minimap2 -y -ax map-ont -k8 -w5 -s 80 --score-N 0 "
            f"--secondary no --sam-hit-only --MD {xfasta_path} - | "
            f"samtools view -F0x800 -bho - | samtools sort -o {sorted_bam}"
        )
        os.system(cmd)
        os.system(f"samtools index {sorted_bam}")
    return sorted_bam

# ---------------------------------------------------------------------
# FASTA & BED helpers
# ---------------------------------------------------------------------
def fasta_to_xfasta(input_fa, ref_dir, prefix="x"):
    if not os.path.isfile(input_fa):
        raise FileNotFoundError(f"FASTA not found: {input_fa}")
    xfasta = os.path.join(ref_dir, prefix + os.path.basename(input_fa))
    os.system(f"python lib/xr_fasta2x_rc.py {input_fa} {xfasta}")
    return xfasta

def sanitize_fasta(input_fa, ref_dir, suffix="_clean"):
    out_fa = os.path.join(ref_dir, Path(input_fa).stem + suffix + ".fa")
    alias_map, counter = {}, 1
    with open(input_fa) as fin, open(out_fa, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                orig = line[1:].strip()
                alias = f"contig{counter}"
                fout.write(f">{alias}\n")
                alias_map[orig] = alias
                counter += 1
            else:
                fout.write(line)
    return out_fa, alias_map

def bed_gen_with_alias(xfasta_file, alias_map, xna_base, sub_base,
                       chunk_range, chunk_shift, output_bed):
    with open(output_bed, "w") as fr, open(xfasta_file) as fo:
        for line in fo:
            if not line.startswith(">"): continue
            header = line[1:].strip()
            alias = alias_map.get(header)
            if not alias: continue
            x_pos_base = fetch_xna_pos(header)
            for xb in x_pos_base:
                x_base = xb[0]
                pos_str = xb[1].replace("-", "").replace("]", "")
                x_pos = int("".join(filter(str.isdigit, pos_str)))
                strand = "+" if x_base == xna_base else "-"
                start = x_pos - chunk_range + chunk_shift
                end   = x_pos + chunk_range + 1 + chunk_shift
                fr.write(f"{alias}\t{start}\t{end}\t{header}\t0\t{strand}\n")
    return output_bed

# ---------------------------------------------------------------------
# Motif extraction (from BED + FASTA)
# ---------------------------------------------------------------------
def extract_context_from_bed(bed_file, fasta_file, flank=4):
    """
    Extract motif centered on XNA from BED + FASTA.
    Example BED: contig1 69 70 ref_PZ_NB25+XPOS[P:69] 0 +
    """
    with open(bed_file) as f:
        line = f.readline().strip()
    if not line:
        raise ValueError(f"Empty BED file: {bed_file}")
    fields = line.split("\t")
    chrom, start, end, header, _, strand = fields
    start, end = int(start), int(end)

    m = re.search(r"\[([A-Za-z]):\d+\]", header)
    if not m:
        raise ValueError(f"Cannot parse XNA from header: {header}")
    xna_base = m.group(1).upper()

    record = next((r for r in SeqIO.parse(fasta_file, "fasta") if r.id == chrom), None)
    if record is None:
        raise ValueError(f"Contig {chrom} not found in {fasta_file}")
    seq = str(record.seq)

    center = start
    motif_start = max(0, center - flank)
    motif_end   = min(len(seq), center + flank + 1)
    motif_seq = seq[motif_start:motif_end]
    center_index = center - motif_start

    motif_seq = motif_seq[:center_index] + xna_base + motif_seq[center_index + 1:]
    print(f"[INFO] Motif extracted: {motif_seq} (center {center_index}) from {chrom}:{start}-{end}")
    return motif_seq, center_index, xna_base

# ---------------------------------------------------------------------
# Basecall-anchored chunk generation
# ---------------------------------------------------------------------

def generate_mod_chunks(pod_file, bam_file, chunk_dir, bed_file, kmer_context,
                        kmer_table_path, fasta_file):
    chunk_file = os.path.join(chunk_dir, "mod_chunks.npz")

    # --- Run chunks only if regenerate_chunks is True ---
    print(f"[DEBUG] regenerate_chunks = {regenerate_chunks}")
    if regenerate_chunks:
        print("Xemora [STATUS] - Generating modified chunks")
    else:
        print("Xemora [STATUS] - Skipping modified chunk generation")
        return chunk_file

    # Handle "4 4" or integer forms gracefully
    try:
        if isinstance(kmer_context, str):
            parts = [int(x) for x in kmer_context.split()]
            flank = parts[0]
        else:
            flank = int(kmer_context) // 2
    except Exception as e:
        raise ValueError(f"Invalid kmer_context format: {kmer_context}") from e

    motif_seq, center_index, xna_base = extract_context_from_bed(bed_file, fasta_file, flank)
    motif_cmds = ' '.join([f"--motif {motif_seq.replace(xna_base, b)} {center_index}" for b in ["A", "T", "G", "C"]])

    cmd = (
        f"remora dataset prepare {pod_file} {bam_file} "
        f"--output-remora-training-file {chunk_file} "
        f"{motif_cmds} "
        f"--mod-base {xna_base} {xna_base} "
        f"--kmer-context-bases {kmer_context} "
        f"--max-chunks-per-read {2*mod_chunk_range+1} "
        f"--refine-kmer-level-table {kmer_table_path} "
        f"--refine-rough-rescale "
        f"--chunk-context {chunk_context} "
        f"--basecall-anchor "
        f"--refine-scale-iters 3 "
        f"--refine-half-bandwidth 80 "
        f"--refine-short-dwell-parameters 8 3 2.0"


    )

    print(f"[DEBUG] Remora (modified):\n{cmd}")
    os.system(cmd)
    return chunk_file



def generate_can_chunks(pod_file, bam_file, chunk_dir, bed_file, kmer_context,
                        kmer_table_path, fasta_file):
    chunk_file = os.path.join(chunk_dir, "can_chunks.npz")

    print(f"[DEBUG] regenerate_chunks = {regenerate_chunks}")
    if regenerate_chunks:
        print("Xemora [STATUS] - Generating canonical chunks")
    else:
        print("Xemora [STATUS] - Skipping canonical chunk generation")
        return chunk_file


    # Handle "4 4" or integer forms gracefully
    try:
        if isinstance(kmer_context, str):
            parts = [int(x) for x in kmer_context.split()]
            flank = parts[0]
        else:
            flank = int(kmer_context) // 2
    except Exception as e:
        raise ValueError(f"Invalid kmer_context format: {kmer_context}") from e

    motif_seq, center_index, xna_base = extract_context_from_bed(bed_file, fasta_file, flank)
    motif_cmds = ' '.join([f"--motif {motif_seq.replace(xna_base, b)} {center_index}" for b in ["A", "T", "G", "C"]])
    print(motif_cmds)

    # Again: all one line, no newline separators
    cmd = (
        f"remora dataset prepare {pod_file} {bam_file} "
        f"--output-remora-training-file {chunk_file} "
        f"{motif_cmds} "
        f"--mod-base-control "
        f"--kmer-context-bases {kmer_context} "
        f"--max-chunks-per-read {2*mod_chunk_range+1} "
        f"--refine-kmer-level-table {kmer_table_path} "
        f"--refine-rough-rescale "
        f"--chunk-context {chunk_context} "
        f"--basecall-anchor "
        f"--refine-scale-iters 3 "
        f"--refine-half-bandwidth 80 "
        f"--refine-short-dwell-parameters 8 3 2.0"

        
    )


    print(f"[DEBUG] Remora (canonical):\n{cmd}")
    os.system(cmd)
    return chunk_file


# ---------------------------------------------------------------------
# Merge + Train
# ---------------------------------------------------------------------
def merge_chunks(chunk_dir, mod_chunks, can_chunks, balance_chunks):
    """Merge modified and canonical Remora datasets for basecall-anchored training."""
    out = os.path.join(chunk_dir, 'training_chunks.npz')
    if not remerge_chunks:
        print('Xemora [STATUS] - Skipping chunk merging')
        return out

    print('Xemora [STATUS] - Merging modified + canonical datasets')

    # Build command (basecall-anchored: take all chunks, not chunk_num_000)
    if balance_chunks:
        cmd = (
            "remora dataset merge "
            f"--balance "
            f"--input-dataset {os.path.join(chunk_dir,'mod_chunks.npz')} all "
            f"--input-dataset {os.path.join(chunk_dir,'can_chunks.npz')} all "
            f"--output-dataset {out}"
        )
    else:
        cmd = (
            "remora dataset merge "
            f"--input-dataset {mod_chunks} all "
            f"--input-dataset {can_chunks} all "
            f"--output-dataset {out}"
        )

    print(f"[DEBUG] Running: {cmd}")
    os.system(cmd)

    # Verify successful merge
    if not os.path.exists(out):
        raise FileNotFoundError(
            f"[ERROR] Merge failed — no file found at {out}. "
            "Check remora output above for errors."
        )

    return out



def xemora_training(model_dir, training_chunks):
    if not gen_model:
        print("Xemora [STATUS] - Skipping training")
        return os.path.join(model_dir, "model_best.pt")
    print("Xemora [STATUS] - Training Xemora model")
    cmd = (
        f"remora model train {training_chunks} "
        f"--model {ml_model_path} "
        f"--device 0 "
        f"--output-path {model_dir} "
        f"--overwrite "
        f"--kmer-context-bases {kmer_context} "
        f"--chunk-context {chunk_context} "
        f"--val-prop {val_proportion} "
        f"--batch-size 100"
    )
    os.system(cmd)
    return os.path.join(model_dir, "model_best.pt")

# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main():
    # Step 0: FASTA + BED prep
    mod_xfasta = fasta_to_xfasta(xna_ref_fasta, ref_dir, "x")
    can_xfasta = fasta_to_xfasta(dna_ref_fasta, ref_dir, "x")
    mod_xfasta_clean, mod_alias = sanitize_fasta(mod_xfasta, ref_dir)
    can_xfasta_clean, can_alias = sanitize_fasta(can_xfasta, ref_dir)

    mod_bed = bed_gen_with_alias(mod_xfasta, mod_alias, mod_base, mod_base,
                                 mod_chunk_range, mod_chunk_shift,
                                 os.path.join(ref_dir, f"{mod_base}.bed"))
    can_bed = bed_gen_with_alias(can_xfasta, can_alias, mod_base, can_base,
                                 can_chunk_range, can_chunk_shift,
                                 os.path.join(ref_dir, f"{can_base}.bed"))

    # Step 1: POD5 prep
    mod_ft = validate_read_directory(xna_raw_dir)
    can_ft = validate_read_directory(dna_raw_dir)
    mod_pod5 = cod5_to_fast5(xna_raw_dir, mod_pod_dir, overwrite_pod) if mod_ft == "fast5" else pod5_merge(xna_raw_dir, mod_pod_dir, overwrite_pod)
    can_pod5 = cod5_to_fast5(dna_raw_dir, can_pod_dir, overwrite_pod) if can_ft == "fast5" else pod5_merge(dna_raw_dir, can_pod_dir, overwrite_pod)

    # Step 2: Basecall
    mod_bam = dorado_basecall(dorado_path, dorado_model, min_qscore, mod_pod5, mod_bam_dir, basecall_pod, max_mod_reads, filter_mod_readIDs)
    can_bam = dorado_basecall(dorado_path, dorado_model, min_qscore, can_pod5, can_bam_dir, basecall_pod, max_can_reads, filter_can_readIDs)

    # Step 3: Align
    mod_aln = minimap2_aligner(mod_bam, mod_xfasta_clean, mod_bam_dir)
    can_aln = minimap2_aligner(can_bam, can_xfasta_clean, can_bam_dir)

    # Step 4: Generate motif-anchored chunks
    mod_chunks = generate_mod_chunks(mod_pod5, mod_aln, chunk_dir, mod_bed, kmer_context, kmer_table_path, mod_xfasta_clean)
    can_chunks = generate_can_chunks(can_pod5, can_aln, chunk_dir, can_bed, kmer_context, kmer_table_path, can_xfasta_clean)

    # Step 5: Merge + Train
    training_chunks = merge_chunks(chunk_dir, mod_chunks, can_chunks, balance_chunks)
    xemora_training(model_dir, training_chunks)

if __name__ == "__main__":
    main()

