########################################################################
########################################################################
"""
xr_train_methods.py

Single-context Xemora model training with Dorado basecalling and Remora.
Normalized via booleans in xr_params:
  - Dataset: USE_PA_SCALING / USE_KMER_REFINE / USE_ROUGH_RESCALE
  - Plots  : PLOT_USE_KMER_REFINE / PLOT_USE_ROUGH_RESCALE

Keeps kmer_table_path relative in xr_params, resolves to absolute at runtime.
"""
########################################################################
########################################################################

import os
import sys
from pathlib import Path
from typing import Optional
from xr_tools  import *          # project helpers (fetch_xna_pos, xna_base_rc, etc.)
from xr_params import *          # project parameters & booleans

print('Xenoquant [STATUS] - Initializing Xenoquant...')

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
# Paths & Dirs
# ---------------------------------------------------------------------
working_dir = check_make_dir(working_dir)
ref_dir     = check_make_dir(os.path.join(working_dir, 'references'))
chunk_dir   = check_make_dir(os.path.join(working_dir, 'chunks'))
model_dir   = check_make_dir(os.path.join(working_dir, 'model'))

mod_dir     = check_make_dir(os.path.join(working_dir, 'modified_preprocess'))
mod_pod_dir = check_make_dir(os.path.join(mod_dir, 'pod5'))
mod_bam_dir = check_make_dir(os.path.join(mod_dir, 'bam'))

can_dir     = check_make_dir(os.path.join(working_dir, 'canonical_preprocess'))
can_pod_dir = check_make_dir(os.path.join(can_dir, 'pod5'))
can_bam_dir = check_make_dir(os.path.join(can_dir, 'bam'))

# Make all relative params resolve relative to this file’s dir, not CWD
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

def resolve_param_path(p: str) -> str:
    """Expand ~ and resolve relative to this file’s directory."""
    return os.path.abspath(os.path.join(BASE_DIR, os.path.expanduser(p)))

# ---------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------
def validate_read_directory(raw_dir):
    """Return 'fast5' or 'pod5' if the directory is homogeneous; else warn."""
    if not os.path.isdir(raw_dir):
        print(f"Xenoquant [ERROR] - Directory not found: {raw_dir}")
        sys.exit(1)
    exts = {f.split('.')[-1] for f in os.listdir(raw_dir) if '.' in f}
    if len(exts) != 1:
        print(f"Xenoquant [STATUS] - Passed reads directory not homogeneous. File types found: {sorted(exts)}")
        return ""
    return list(exts)[0]

def cod5_to_fast5(fast5_input, pod_dir, overwrite_pod):
    """Convert directory of FAST5 to a single POD5."""
    out = os.path.join(pod_dir, os.path.basename(fast5_input) + '.pod5')
    if overwrite_pod or not os.path.exists(out):
        cmd = f'pod5 convert fast5 --force-overwrite {fast5_input}/*.fast5 -o {out}'
        os.system(cmd)
    else:
        print('Xenoquant [STATUS] - Skipping FAST5→POD5 conversion')
    return out

def pod5_merge(pod5_input, pod_dir, overwrite_pod):
    """Merge many POD5 into one."""
    out = os.path.join(pod_dir, os.path.basename(pod5_input) + '.pod5')
    if overwrite_pod or not os.path.exists(out):
        cmd = f'pod5 merge --force-overwrite {pod5_input}/*.pod5 -o {out}'
        os.system(cmd)
    else:
        print('Xenoquant [STATUS] - Skipping POD5 merge')
    return out

def dorado_basecall(dorado_path, dorado_model, min_qscore, pod5_file, bam_dir, basecall_pod, max_reads, filter_readIDs):
    """Run Dorado basecaller → BAM."""
    out_bam = os.path.join(bam_dir, 'bc.bam')
    if basecall_pod or not os.path.exists(out_bam):
        print('Xenoquant [STATUS] - Performing basecalling using Dorado')
        args = f'--no-trim --emit-moves {pod5_file}'
        if filter_readIDs: args += f' -l {filter_readIDs}'
        if max_reads:      args += f' -n {max_reads}'
        if min_qscore:     args += f' --min-qscore {min_qscore}'
        cmd = f'{os.path.expanduser(dorado_path)} basecaller {dorado_model} {args} > {out_bam}'
        os.system(cmd)
    else:
        print('Xenoquant [STATUS] - Skipping basecalling; BAM exists.')
    return out_bam

def minimap2_aligner(input_bam, xfasta_path, bam_directory):
    """Align with minimap2 → sorted + indexed BAM."""
    raw_bam   = os.path.join(bam_directory, "aligned.unsorted.bam")
    sorted_bam= os.path.join(bam_directory, "aligned.sorted.bam")
    if regenerate_bam:
        cmd_align = (
            f"samtools fastq -T '*' {input_bam} | "
            f"minimap2 -y -ax map-ont -k8 -w5 -s 80 --score-N 0 "
            f"--secondary no --sam-hit-only --MD {xfasta_path} - | "
            f"samtools view -F0x800 -bho {raw_bam}"
        )
        os.system(cmd_align)
        os.system(f"samtools sort -o {sorted_bam} {raw_bam}")
        os.system(f"samtools index {sorted_bam}")
        if os.path.exists(raw_bam):
            os.remove(raw_bam)
            print(f"[CLEANUP] Removed unsorted BAM: {raw_bam}")
    else:
        print('Xenoquant [STATUS] - Skipping Minimap2 alignment.')
    return sorted_bam

def ensure_index(bam_path):
    if not (os.path.exists(bam_path + ".bai") or os.path.exists(bam_path + ".csi")):
        print(f"[INFO] Index missing, indexing {bam_path}")
        os.system(f"samtools index {bam_path}")
    else:
        print(f"[INFO] Index present for {bam_path}")

# ---------------------------------------------------------------------
# FASTA/xFASTA/BED handling
# ---------------------------------------------------------------------
def fasta_to_xfasta(input_fa, ref_dir, prefix="x"):
    """Convert FASTA → xFASTA via project converter."""
    if not os.path.isfile(input_fa):
        raise FileNotFoundError(f"[ERROR] FASTA not found: {input_fa}")
    xfasta = os.path.join(ref_dir, prefix + os.path.basename(input_fa))
    cmd = f'python lib/xr_fasta2x_rc.py {os.path.expanduser(input_fa)} {xfasta}'
    rc = os.system(cmd)
    if rc != 0 or not os.path.exists(xfasta):
        raise RuntimeError(f"[ERROR] xFASTA conversion failed: {cmd}")
    return xfasta

def rc_xfasta_candidate(path: str) -> Optional[str]:
    """Return candidate RC filename our converter writes, if it exists."""
    if path.endswith(".fasta"):
        cand = path[:-6] + "_rc.fasta"
    elif path.endswith(".fa"):
        cand = path[:-3] + "_rc.fa"
    else:
        cand = path + "_rc.fa"
    return cand if os.path.exists(cand) else None

def xfasta_requires_rc(xfasta_path: str, xna_base: str) -> bool:
    """Return True if headers indicate XNA on reverse strand."""
    rc_symbol = xna_base_rc(xna_base, xna_base_pairs)
    with open(os.path.expanduser(xfasta_path)) as f:
        for ln in f:
            if ln.startswith(">"):
                for xb in fetch_xna_pos(ln[1:].strip()):
                    base = xb[0]
                    if base == rc_symbol:
                        return True
    return False

def sanitize_fasta(input_fa, ref_dir, suffix="_clean"):
    """Rewrite headers to BAM-safe aliases: contig1, contig2, ..."""
    out_fa = os.path.join(ref_dir, Path(input_fa).stem + suffix + ".fa")
    map_tsv= os.path.join(ref_dir, "contig_map.tsv")
    alias_map, counter = {}, 1
    with open(input_fa) as fin, open(out_fa, "w") as fout, open(map_tsv, "w") as fmap:
        for line in fin:
            if line.startswith(">"):
                orig = line[1:].strip()
                alias = f"contig{counter}"; counter += 1
                alias_map[orig] = alias
                fout.write(f">{alias}\n")
                fmap.write(f"{alias}\t{orig}\n")
            else:
                fout.write(line)
    print(f"Xenoquant [STATUS] - Wrote sanitized FASTA: {out_fa}")
    print(f"Xenoquant [STATUS] - Wrote mapping TSV: {map_tsv}")
    return out_fa, alias_map

def bed_gen_with_alias(xfasta_file, alias_map, xna_base, sub_base,
                       chunk_range, chunk_shift, output_bed):
    """Generate BED with alias in col1 and original header in col4."""
    with open(output_bed, "w") as fr, open(xfasta_file) as fo:
        for line in fo:
            if not line.startswith(">"): continue
            header = line[1:].strip()
            alias  = alias_map.get(header)
            if not alias:
                print(f"[WARN] No alias for header: {header}")
                continue
            x_pos_base = fetch_xna_pos(header)
            if not x_pos_base:
                print(f"[WARN] No XNA positions in header: {header}")
                continue
            for xb in x_pos_base:
                x_base, pos_str = xb[0], (xb[1] if len(xb) > 1 else xb[1:])
                try:
                    x_pos = int("".join(filter(str.isdigit, pos_str)))
                except Exception:
                    print(f"[WARN] Could not parse position from {pos_str} in {header}")
                    continue
                strand = "+" if x_base == xna_base else "-"
                start  = x_pos - chunk_range + chunk_shift
                end    = x_pos + chunk_range + 1 + chunk_shift
                fr.write(f"{alias}\t{start}\t{end}\t{header}\t0\t{strand}\n")
    print(f"Xenoquant [STATUS] - Wrote BED: {output_bed}")
    return output_bed





# ---------------------------------------------------------------------
# Chunk prep / merge / train
# ---------------------------------------------------------------------
def generate_mod_chunks(pod_file, bam_file, chunk_dir, bed_file, mod_base, kmer_context, kmer_table_path, regenerate_chunks):
    """
    Prepares and runs a Remora command to generate chunks for modified base analysis.

    Parameters:
    pod_file (str): Merged POD5 file path.
    bam_file (str): BAM file generated from POD5 file basecalling.
    bed_file (str): BED file specifying the regions to analyze.
    mod_base (str): Modified base to perform inference on, defined in xr_params.
    kmer_context (str): k-mer context to use, defined in xr_params.
    kmer_table_path (str): Path to the k-mer table, defined in xr_params.
    regenerate_chunks (bool): Parameter to allow chunk files to be regenerated.

    Returns: 
    str: File path to the generated Remora chunk file.
    """
    chunk_file = os.path.join(chunk_dir, 'mod_chunks.npz')

    if regenerate_chunks:
        print('Xenoquant [STATUS] - Generating chunks for modified basecalling.')

        cmd = (
            f"remora dataset prepare "
            f"{pod_file} "
            f"{bam_file} "
            f"--output-remora-training-file {chunk_file} "
            f"--focus-reference-positions {bed_file} "
            f"--mod-base {mod_base} {mod_base} "
            f"--kmer-context-bases {kmer_context} "
            f"--max-chunks-per-read {2 * mod_chunk_range + 1} "
            f"--refine-kmer-level-table {kmer_table_path} "
            f"--refine-rough-rescale "
            f"--refine-scale-iters -1 "
            f"--motif {can_base} 0 "
            f"--chunk-context {chunk_context}"
        )

        os.system(cmd)
        return chunk_file

    else:
        print('Xenoquant [STATUS] - Skipping modified chunk generation')
        return chunk_file
        
'''
            f"--refine-half-bandwidth 80 "
            f"--refine-short-dwell-parameters 8 3 2.0 "
'''


def generate_can_chunks(pod_file, bam_file, chunk_dir, bed_file, can_base, kmer_context, kmer_table_path, regenerate_chunks):
    """
    Prepares and runs a Remora command to generate chunks for canonical base analysis.

    Parameters:
    pod_file (str): Merged POD5 file path.
    bam_file (str): BAM file generated from POD5 file basecalling.
    bed_file (str): BED file specifying the regions to analyze.
    can_base (str): Canonical base to perform inference on, defined in xr_params.
    kmer_context (str): k-mer context to use, defined in xr_params.
    kmer_table_path (str): Path to the k-mer table, defined in xr_params.
    regenerate_chunks (bool): Parameter to allow chunk files to be regenerated.

    Returns: 
    str: File path to the generated Remora chunk file.
    """
    chunk_file = os.path.join(chunk_dir, 'can_chunks.npz')

    if regenerate_chunks:
        print('Xenoquant [STATUS] - Generating chunks for modified basecalling.')

        cmd = (
            f"remora dataset prepare "
            f"{pod_file} "
            f"{bam_file} "
            f"--output-remora-training-file {chunk_file} "
            f"--focus-reference-positions {bed_file} "
            f"--mod-base-control "
            f"--kmer-context-bases {kmer_context} "
            f"--max-chunks-per-read {2 * mod_chunk_range + 1} "
            f"--refine-kmer-level-table {kmer_table_path} "
            f"--refine-rough-rescale "
            f"--refine-scale-iters -1 "
            f"--motif {can_base} 0 "
            f"--chunk-context {chunk_context}"
        )


        os.system(cmd)

    else:
        print('Xenoquant [STATUS] - Skipping canonical chunk generation')
    
    return chunk_file


def merge_chunks(chunk_dir, mod_chunks, can_chunks, balance_chunks):
    out = os.path.join(chunk_dir, 'training_chunks.npz')
    if not remerge_chunks:
        print('Xenoquant [STATUS] - Skipping chunk merging'); return out
    print('Xenoquant [STATUS] - Merging chunks for training.')
    if balance_chunks:
        cmd = (
            "remora dataset merge "
            f"--balance "
            f"--input-dataset {os.path.join(chunk_dir,'mod_chunks.npz')} {chunk_num}_000 "
            f"--input-dataset {os.path.join(chunk_dir,'can_chunks.npz')} {chunk_num}_000 "
            f"--output-dataset {out}"
        )
    else:
        cmd = (
            "remora dataset merge "
            f"--input-dataset {mod_chunks} {chunk_num}_000 "
            f"--input-dataset {can_chunks} {chunk_num}_000 "
            f"--output-dataset {out}"
        )
    os.system(cmd)
    return out

def Xenoquant_training(model_dir, training_chunks):
    if not gen_model:
        print('Xenoquant [STATUS] - Skipping model training')
        return os.path.join(model_dir, 'model_best.pt'), os.path.join(model_dir, 'validation.log')
    print('Xenoquant [STATUS] - Training model.')
    cmd = (
        "remora model train "
        f"{os.path.join(chunk_dir,'training_chunks.npz')} "
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
    model_path = os.path.join(model_dir, 'model_best.pt')
    val_log    = os.path.join(model_dir, 'validation.log')
    return model_path, val_log

import subprocess
import shlex
import os

def run_remora_ref_region_plot(
    can_pod5, can_bam,
    mod_pod5, mod_bam,
    ref_bed,
    highlight_bed,
    levels_table,
    out_dir,
    prefix="ref_region",
    log_name="ref_region.log",
    soft_fail=True,                 # <— NEW
):
    """
    Generate Remora ref_region plots. If soft_fail=True, any error is logged and
    the pipeline continues (returns a nonzero rc). If False, raises on error.
    """
    os.makedirs(out_dir, exist_ok=True)
    # Resolve absolute paths (keeps existing behavior)
    ref_bed      = os.path.abspath(os.path.expanduser(ref_bed))
    highlight_bed= os.path.abspath(os.path.expanduser(highlight_bed))
    levels_table = os.path.abspath(os.path.expanduser(levels_table))

    if not os.path.exists(levels_table):
        msg = f"[ERROR] k-mer levels table not found: {levels_table}"
        if soft_fail:
            print(msg)
            return 2  # arbitrary nonzero
        raise FileNotFoundError(msg)

    cmd = (
        "remora analyze plot ref_region "
        f"--pod5-and-bam {can_pod5} {can_bam} "
        f"--pod5-and-bam {mod_pod5} {mod_bam} "
        f"--ref-regions {ref_bed} "
        f"--highlight-ranges {highlight_bed} "
        f"--refine-kmer-level-table {levels_table} "
        f"--refine-rough-rescale "
        f"--log-filename {log_name}"
    )
    print(f"[DEBUG] Running Remora ref_region plot:\n{cmd}")

    try:
        rc = subprocess.run(shlex.split(cmd), cwd=out_dir).returncode
    except Exception as e:
        if soft_fail:
            print(f"[WARN] Remora ref_region plotting raised {type(e).__name__}: {e}")
            return 3
        raise

    if rc != 0:
        msg = f"[WARN] Remora ref_region plotting failed with code {rc}. " \
              f"Continuing (soft_fail=True). See {os.path.join(out_dir, log_name)}"
        if soft_fail:
            print(msg)
            return rc
        raise RuntimeError(msg)

    print(f"Xenoquant [STATUS] - Remora ref_region plots written to {out_dir}")
    return 0



# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main():
    # Step 0: xFASTA conversion
    mod_xfasta = fasta_to_xfasta(xna_ref_fasta, ref_dir, prefix="x")
    can_xfasta = fasta_to_xfasta(dna_ref_fasta, ref_dir, prefix="x")

    # Strand/orientation
    if xfasta_requires_rc(mod_xfasta, mod_base):
        cand_mod = rc_xfasta_candidate(mod_xfasta)
        cand_can = rc_xfasta_candidate(can_xfasta)
        if not cand_mod or not cand_can:
            raise FileNotFoundError(f"[ERROR] RC xFASTA expected but not found for mod:{mod_xfasta}, can:{can_xfasta}")
        print(f"[INFO] XNA on reverse strand; using RC xFASTAs: {cand_mod}, {cand_can}")
        mod_xfasta, can_xfasta = cand_mod, cand_can
    else:
        print("[INFO] XNA on forward strand; using forward xFASTAs")

    # Sanitize after deciding on RC
    mod_xfasta_clean, mod_alias_map = sanitize_fasta(mod_xfasta, ref_dir)
    can_xfasta_clean, can_alias_map = sanitize_fasta(can_xfasta, ref_dir)

    # BEDs from the SAME FASTAs used for alignment
    mod_bed_file = bed_gen_with_alias(
        mod_xfasta, mod_alias_map, mod_base, mod_base,
        mod_chunk_range, mod_chunk_shift,
        os.path.join(ref_dir, f"{mod_base}.bed")
    )
    can_bed_file = bed_gen_with_alias(
        can_xfasta, can_alias_map, mod_base, can_base,
        can_chunk_range, can_chunk_shift,
        os.path.join(ref_dir, f"{can_base}.bed")
    )
    # Ref regions (for plots; use whichever you prefer)
    ref_regions_bed = bed_gen_with_alias(
        mod_xfasta, mod_alias_map, mod_base, mod_base,
        FLANK, 0,
        os.path.join(ref_dir, f"{mod_base}_ref_regions.bed")
    )

    # Step 1: POD5 prep
    mod_ft = validate_read_directory(xna_raw_dir)
    can_ft = validate_read_directory(dna_raw_dir)
    if mod_ft not in ("fast5", "pod5") or can_ft not in ("fast5", "pod5"):
        print('Xenoquant [ERROR] - Raw data directory must contain only POD5 or FAST5'); sys.exit(1)
    mod_merged_pod5 = cod5_to_fast5(xna_raw_dir, mod_pod_dir, overwrite_pod) if mod_ft == 'fast5' else pod5_merge(xna_raw_dir, mod_pod_dir, overwrite_pod)
    can_merged_pod5 = cod5_to_fast5(dna_raw_dir, can_pod_dir, overwrite_pod) if can_ft == 'fast5' else pod5_merge(dna_raw_dir, can_pod_dir, overwrite_pod)

    # Step 2: Basecall
    mod_bc_bam = dorado_basecall(dorado_path, dorado_model, min_qscore, mod_merged_pod5, mod_bam_dir, basecall_pod, max_mod_reads, filter_mod_readIDs)
    can_bc_bam = dorado_basecall(dorado_path, dorado_model, min_qscore, can_merged_pod5, can_bam_dir, basecall_pod, max_can_reads, filter_can_readIDs)

    # Step 3: Align against SANITIZED xFASTAs
    mod_aligned_bam = minimap2_aligner(mod_bc_bam, mod_xfasta_clean, mod_bam_dir)
    can_aligned_bam = minimap2_aligner(can_bc_bam, can_xfasta_clean, can_bam_dir)
    ensure_index(mod_aligned_bam)
    ensure_index(can_aligned_bam)

    # -----------------------------------------------------------------
    # Optional: Generate Remora ref_region plots
    # -----------------------------------------------------------------
    if generate_remora_plots:
        remora_plot_dir = check_make_dir(os.path.join(working_dir, "remora_plots"))
        rc_plot = run_remora_ref_region_plot(
            can_merged_pod5, can_aligned_bam,
            mod_merged_pod5, mod_aligned_bam,
            ref_bed=ref_regions_bed,
            highlight_bed=mod_bed_file,
            levels_table=kmer_table_path,
            out_dir=remora_plot_dir,
            soft_fail=True,   # continue even if plotting fails
        )
        if rc_plot != 0:
            print(f"[INFO] Skipping plot-dependent post-steps; plot rc={rc_plot}")
    else:
        print("Xenoquant [STATUS] - Skipping Remora plotting (generate_remora_plots=False)")




    # Step 4: Generate Remora chunks
    mod_chunk_path = generate_mod_chunks(mod_merged_pod5, mod_aligned_bam, chunk_dir, mod_bed_file, mod_base, kmer_context, kmer_table_path, regenerate_chunks)
    can_chunk_path = generate_can_chunks(can_merged_pod5, can_aligned_bam, chunk_dir, can_bed_file, can_base, kmer_context, kmer_table_path, regenerate_chunks)

    # Step 5: Merge
    training_chunk_path = merge_chunks(chunk_dir, mod_chunk_path, can_chunk_path, balance_chunks)

    # Step 6: Train
    model_path, validation_log_path = Xenoquant_training(model_dir, training_chunk_path)

if __name__ == "__main__":
    main()

