#!/usr/bin/env python3
########################################################################
"""
xr_train_methods.py

Single-context Xemora model training pipeline with correct order:
FASTA -> xFASTA (annotated) -> BED -> sanitize xFASTA/BED -> alignment -> plotting -> chunks -> train
"""
########################################################################

import os
import sys
from pathlib import Path
from xr_tools  import *
from xr_params import *

print('Xemora [Status] - Initializing Xemora...')

########################################################################
# Args & Dirs
########################################################################

working_dir   = os.path.expanduser(sys.argv[1])
xna_raw_dir   = os.path.expanduser(sys.argv[2])
xna_ref_fasta = os.path.expanduser(sys.argv[3])
dna_raw_dir   = os.path.expanduser(sys.argv[4])
dna_ref_fasta = os.path.expanduser(sys.argv[5])

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

########################################################################
# Utilities
########################################################################

def validate_read_directory(raw_dir):
    """Check read dir is homogeneous FAST5 or POD5."""
    if not os.path.isdir(raw_dir):
        print(f'Xemora [ERROR] - Input reads directory does not exist: {raw_dir}')
        sys.exit(1)
    exts = {p.split('.')[-1] for p in os.listdir(raw_dir) if not p.startswith('.')}
    if len(exts) != 1:
        print(f'Xemora [STATUS] - Directory not homogeneous. File types: {list(exts)}')
        return ""
    return list(exts)[0]

def cod5_to_fast5(fast5_input, pod_dir, overwrite_pod):
    """Convert FAST5 → POD5."""
    pod5_output = os.path.join(pod_dir, os.path.basename(fast5_input) + '.pod5')
    if overwrite_pod or not os.path.exists(pod5_output):
        cmd = f'pod5 convert fast5 --force-overwrite {fast5_input}/*.fast5 -o {pod5_output}'
        os.system(cmd)
    else:
        print('Xemora [STATUS] - Skipping FAST5 conversion')
    return pod5_output

def pod5_merge(pod5_input, pod_dir, overwrite_pod):
    """Merge POD5 files."""
    merged_pod5 = os.path.join(pod_dir, os.path.basename(pod5_input) + '.pod5')
    if overwrite_pod or not os.path.exists(merged_pod5):
        cmd = f'pod5 merge --force-overwrite {pod5_input}/*.pod5 -o {merged_pod5}'
        os.system(cmd)
    else:
        print('Xemora [STATUS] - Skipping POD5 merge')
    return merged_pod5

def dorado_basecall(dorado_path, dorado_model, min_qscore, pod_input, bam_directory,
                    basecall_pod, max_reads, filter_readIDs):
    """Basecall with Dorado (emit moves, no trim)."""
    output_bam = os.path.join(bam_directory, 'bc.bam')
    if basecall_pod or not os.path.exists(output_bam):
        print('Xemora [STATUS] Performing basecalling with Dorado')
        dorado_bin = os.path.expanduser(dorado_path)
        dorado_args = f'--no-trim --emit-moves {pod_input}'
        if filter_readIDs: dorado_args += f' -l {filter_readIDs}'
        if max_reads:      dorado_args += f' -n {max_reads}'
        if min_qscore:     dorado_args += f' --min-qscore {min_qscore}'
        cmd = f'{dorado_bin} basecaller {dorado_model} {dorado_args} > {output_bam}'
        os.system(cmd)
    else:
        print('Xemora [STATUS] - Skipping basecalling')
    return output_bam

def minimap2_aligner(input_bam, xfasta_path, bam_directory):
    """Align with minimap2 → sorted BAM."""
    raw_aligned_bam = os.path.join(bam_directory, "aligned.unsorted.bam")
    sorted_bam = os.path.join(bam_directory, "aligned.sorted.bam")
    cmd_align = (
        f"samtools fastq -T '*' {input_bam} | "
        f"minimap2 -y -ax map-ont -k8 -w5 -s 80 --score-N 0 "
        f"--secondary no --sam-hit-only --MD {xfasta_path} - | "
        f"samtools view -F0x800 -bho {raw_aligned_bam}"
    )
    os.system(cmd_align)
    os.system(f"samtools sort -o {sorted_bam} {raw_aligned_bam}")
    os.system(f"samtools index {sorted_bam}")
    if os.path.exists(raw_aligned_bam):
        os.remove(raw_aligned_bam)
        print(f"[CLEANUP] Removed unsorted BAM: {raw_aligned_bam}")
    return sorted_bam



########################################################################
# FASTA + BED Handling
########################################################################

def fasta_to_xfasta(input_fa, ref_dir, prefix="x"):
    """
    Convert a FASTA into xFASTA using the project converter that writes headers
    with XNA focus info (required for fetch_xna_pos).
    """
    if not os.path.isfile(input_fa):
        raise FileNotFoundError(f"[ERROR] FASTA not found: {input_fa}")
    xfasta = os.path.join(ref_dir, prefix + os.path.basename(input_fa))
    cmd = f'python lib/xr_fasta2x_rc.py {os.path.expanduser(input_fa)} {xfasta}'
    rc = os.system(cmd)
    if rc != 0 or not os.path.exists(xfasta):
        raise RuntimeError(f"[ERROR] xFASTA conversion failed: {cmd}")
    return xfasta
    
def rc_xfasta_candidate(path: str) -> str | None:
    """Return the rc xFASTA filename your converter writes, if it exists."""
    if path.endswith(".fasta"):
        cand = path[:-6] + "_rc.fasta"
    elif path.endswith(".fa"):
        cand = path[:-3] + "_rc.fa"
    else:
        cand = path + "_rc.fa"
    return cand if os.path.exists(cand) else None

def xfasta_requires_rc(xfasta_path: str, xna_base: str) -> bool:
    """
    Scan xFASTA headers; if any XPOS base equals xna_base_rc, return True.
    Uses your fetch_xna_pos and xna_base_rc from xr_tools/xr_params.
    """
    rc_symbol = xna_base_rc(xna_base, xna_base_pairs)  # from your params
    with open(os.path.expanduser(xfasta_path)) as f:
        for ln in f:
            if ln.startswith(">"):
                x_pos_base = fetch_xna_pos(ln[1:].strip())
                for xb in x_pos_base:
                    base = xb[0] if isinstance(xb, (list, tuple)) else xb[0]
                    if base == rc_symbol:
                        return True
    return False





def sanitize_fasta(input_fa, ref_dir, suffix="_clean"):
    """
    Rewrite FASTA headers to BAM-safe aliases: contig1, contig2, ...
    Writes:
      - sanitized FASTA (with contigN headers)
      - mapping TSV (alias → original header)
    """
    out_fa = os.path.join(ref_dir, Path(input_fa).stem + suffix + ".fa")
    map_tsv = os.path.join(ref_dir, "contig_map.tsv")

    alias_map = {}
    counter = 1

    with open(input_fa) as fin, open(out_fa, "w") as fout, open(map_tsv, "w") as fmap:
        for line in fin:
            if line.startswith(">"):
                orig = line[1:].strip()
                alias = f"contig{counter}"
                alias_map[orig] = alias
                counter += 1
                fout.write(f">{alias}\n")
                fmap.write(f"{alias}\t{orig}\n")
            else:
                fout.write(line)

    print(f"Xemora [STATUS] - Wrote sanitized FASTA: {out_fa}")
    print(f"Xemora [STATUS] - Wrote mapping TSV: {map_tsv}")
    return out_fa, alias_map


def bed_gen_with_alias(xfasta_file, alias_map, xna_base, sub_base,
                       chunk_range, chunk_shift, output_bed):
    """
    Generate BED file:
      - col1 = BAM-safe alias
      - col4 = original descriptive header
    """
    with open(output_bed, "w") as fr, open(xfasta_file) as fo:
        for line in fo:
            if line.startswith(">"):
                header = line[1:].strip()
                alias = alias_map.get(header)
                if alias is None:
                    print(f"[WARNING] No alias found for header: {header}")
                    continue

                x_pos_base = fetch_xna_pos(header)
                if not x_pos_base:
                    print(f"[WARNING] No XNA positions found in header: {header}")
                    continue

                for x in x_pos_base:
                    if isinstance(x, (tuple, list)) and len(x) >= 2:
                        x_base, pos_str = x[0], x[1]
                    else:
                        x_base, pos_str = x[0], x[1:]  # fallback

                    try:
                        x_pos = int("".join(filter(str.isdigit, pos_str)))
                    except ValueError:
                        print(f"[WARNING] Could not parse position from {pos_str} in {header}")
                        continue

                    strand = "+" if x_base == xna_base else "-"
                    start  = x_pos - chunk_range + chunk_shift
                    end    = x_pos + chunk_range + 1 + chunk_shift
                    fr.write(f"{alias}\t{start}\t{end}\t{header}\t0\t{strand}\n")

    print(f"Xemora [STATUS] - Wrote BED: {output_bed}")
    return output_bed


########################################################################
# Remora wrappers
########################################################################

def plot_ref_region_pair(
    can_pod5, can_bam,
    mod_pod5, mod_bam,
    ref_regions_bed,
    highlight_bed,
    levels_table,
    out_dir,
    log_name="ref_region.log"
):
    """
    remora analyze plot ref_region (rescaling ON).
    No --output-dir; we chdir into out_dir. Use absolute paths for inputs.
    """
    os.makedirs(out_dir, exist_ok=True)
    cwd = os.getcwd()
    try:
        ref_regions_bed = os.path.abspath(os.path.expanduser(ref_regions_bed))
        highlight_bed   = os.path.abspath(os.path.expanduser(highlight_bed))
        levels_table    = os.path.abspath(os.path.expanduser(levels_table))
        if not os.path.exists(levels_table):
            raise FileNotFoundError(f"[ERROR] k-mer levels table not found: {levels_table}")

        os.chdir(out_dir)
        cmd = (
            "remora analyze plot ref_region "
            f"--pod5-and-bam {can_pod5} {can_bam} "
            f"--pod5-and-bam {mod_pod5} {mod_bam} "
            f"--ref-regions {ref_regions_bed} "
            f"--highlight-ranges {highlight_bed} "
            f"--refine-kmer-level-table {levels_table} "
            f"--refine-rough-rescale "
            f"--log-filename {log_name}"
        )
        print(f"[DEBUG] Running Remora ref_region plot:\n{cmd}")
        rc = os.system(cmd)
        if rc != 0:
            raise RuntimeError(f"Remora ref_region plotting failed with code {rc}")
        print(f"Xemora [STATUS] - ref_region plots written to {out_dir}")
    finally:
        os.chdir(cwd)

def prepare_chunks(pod_file, bam_file, chunk_path, can_base,
                   mod_symbol=None, kmer_context=None,
                   kmer_table_path="", is_control=False,
                   max_chunks_per_read=None):
    """
    remora dataset prepare with rescaling always enabled.
    Resolve kmer table to absolute to avoid CWD issues.
    """
    if can_base not in {"A", "C", "G", "T"}:
        raise ValueError(f"can_base must be A/C/G/T, got '{can_base}'")

    kmer_table_abs = os.path.abspath(os.path.expanduser(kmer_table_path))
    if not os.path.exists(kmer_table_abs):
        raise FileNotFoundError(f"k-mer level table not found: {kmer_table_abs}")

    max_chunks = f"--max-chunks-per-read {max_chunks_per_read} " if max_chunks_per_read else ""
    base_cmd = (
        f"remora dataset prepare {pod_file} {bam_file} "
        f"--output-path {chunk_path} --overwrite "
        f"--motif {can_base} 0 "
        f"--kmer-context-bases {kmer_context} "
        f"{max_chunks}"
        f"--refine-kmer-level-table {kmer_table_abs} --refine-rough-rescale "
        f"--chunk-context {chunk_context}"
    )
    if is_control:
        cmd = base_cmd + " --mod-base-control"
    else:
        if not mod_symbol:
            raise ValueError("mod_symbol must be provided for modified dataset.")
        cmd = base_cmd + f" --mod-base {mod_symbol} {can_base}"

    print(f"Xemora [STATUS] - Preparing chunks -> {chunk_path}")
    os.system(cmd)
    return chunk_path

def make_train_config(can_chunks, mod_chunks, chunk_dir, weights=(1,1), name="train_dataset.jsn"):
    """Compose datasets into config for training."""
    cfg_path = os.path.join(chunk_dir, name)
    log_path = os.path.join(chunk_dir, "train_dataset.log")
    w1, w2 = weights
    cmd = (
        f"remora dataset make_config {cfg_path} {can_chunks} {mod_chunks} "
        f"--dataset-weights {w1} {w2} "
        f"--log-filename {log_path}"
    )
    print("Xemora [STATUS] - Creating train dataset config")
    os.system(cmd)
    return cfg_path

def train_model_from_config(config_path, model_dir, model_script,
                            batch_size=512, device=0, val_prop=0.2,
                            lr=0.001, epochs=50, max_updates_per_epoch=10000):
    """Train Remora model with flexible options."""
    os.makedirs(model_dir, exist_ok=True)
    cmd = (
        f"remora model train {config_path} "
        f"--model {model_script} "
        f"--device {device} "
        f"--chunk-context {chunk_context} "
        f"--output-path {model_dir} "
        f"--batch-size {batch_size} "
        f"--lr {lr} "
        f"--epochs {epochs} "
        f"--overwrite "
    )

    print("Xemora [STATUS] - Training model from config")
    print(f"[DEBUG] Running: {cmd}")
    os.system(cmd)
    return os.path.join(model_dir, "model_best.pt"), os.path.join(model_dir, "validation.log")


########################################################################
# Main
########################################################################

def main():
        # Step 0: xFASTA conversion
    mod_xfasta = fasta_to_xfasta(xna_ref_fasta, ref_dir, prefix="x")
    can_xfasta = fasta_to_xfasta(dna_ref_fasta, ref_dir, prefix="x")

    # Swap to RC if required
    if xfasta_requires_rc(mod_xfasta, mod_base):
        cand_mod = rc_xfasta_candidate(mod_xfasta)
        cand_can = rc_xfasta_candidate(can_xfasta)
        if not cand_mod or not cand_can:
            raise FileNotFoundError(
                f"[ERROR] RC xFASTA expected but not found for mod:{mod_xfasta}, can:{can_xfasta}"
            )
        print(f"[INFO] XNA on reverse strand; using RC xFASTAs: {cand_mod}, {cand_can}")
        mod_xfasta = cand_mod
        can_xfasta = cand_can
    else:
        print("[INFO] XNA on forward strand; using forward xFASTAs")


    # Step 0.5: Sanitize *after deciding on RC*
    mod_xfasta_clean, mod_alias_map = sanitize_fasta(mod_xfasta, ref_dir)
    can_xfasta_clean, can_alias_map = sanitize_fasta(can_xfasta, ref_dir)


    # Step 1: BEDs must come from the same FASTAs used for alignment
    mod_bed_file = bed_gen_with_alias(
        mod_xfasta, mod_alias_map,
        mod_base, mod_base,
        mod_chunk_range, mod_chunk_shift,
        os.path.join(ref_dir, f"{mod_base}.bed")
    )

    can_bed_file = bed_gen_with_alias(
        can_xfasta, can_alias_map,
        mod_base, can_base,
        can_chunk_range, can_chunk_shift,
        os.path.join(ref_dir, f"{can_base}.bed")
    )

    # Step 1.5: Ref regions from same RC/forward FASTA
    ref_regions_bed = bed_gen_with_alias(
        mod_xfasta, mod_alias_map,
        mod_base, mod_base,
        flank_size, 0,
        os.path.join(ref_dir, f"{mod_base}_ref_regions.bed")
    )



    # 3) POD5 merge / basecalling (unchanged)
    mod_ft = validate_read_directory(xna_raw_dir)
    can_ft = validate_read_directory(dna_raw_dir)
    if mod_ft not in ("fast5", "pod5") or can_ft not in ("fast5", "pod5"):
        print('Xemora [ERROR] - Raw data directory must contain only POD5 or FAST5'); sys.exit(1)
    mod_merged_pod5 = cod5_to_fast5(xna_raw_dir, mod_pod_dir, overwrite_pod) if mod_ft == 'fast5' else pod5_merge(xna_raw_dir, mod_pod_dir, overwrite_pod)
    can_merged_pod5 = cod5_to_fast5(dna_raw_dir, can_pod_dir, overwrite_pod) if can_ft == 'fast5' else pod5_merge(dna_raw_dir, can_pod_dir, overwrite_pod)

    mod_bc_bam = dorado_basecall(dorado_path, dorado_model, min_qscore, mod_merged_pod5, mod_bam_dir, basecall_pod, max_mod_reads, filter_mod_readIDs)
    can_bc_bam = dorado_basecall(dorado_path, dorado_model, min_qscore, can_merged_pod5, can_bam_dir, basecall_pod, max_can_reads, filter_can_readIDs)

    # 4) Align against SANITIZED xFASTAs
    mod_aligned_bam = minimap2_aligner(mod_bc_bam, mod_xfasta_clean, mod_bam_dir)
    can_aligned_bam = minimap2_aligner(can_bc_bam, can_xfasta_clean, can_bam_dir)


    # 5) Plot (use sanitized BEDs), always rescale (kmer table resolved inside)
    sig_plot_dir = check_make_dir(os.path.join(working_dir, "sig_plots"))
    plot_ref_region_pair(
        can_merged_pod5, can_aligned_bam,
        mod_merged_pod5, mod_aligned_bam,
        ref_regions_bed=ref_regions_bed,
        highlight_bed=can_bed_file,   # << was can_bed (undefined); use canonical BED
        levels_table=kmer_table_path,
        out_dir=sig_plot_dir,
        log_name="ref_region.log"
    )



    # 6) Chunks
    mod_chunk_path = os.path.join(chunk_dir, 'mod_chunks.npz')
    can_chunk_path = os.path.join(chunk_dir, 'can_chunks.npz')
    mod_symbol = mod_label if 'mod_label' in globals() else mod_base
    if can_base not in {"A","C","G","T"}:
        raise ValueError(f"can_base must be A/C/G/T, got {can_base}")

    prepare_chunks(
        mod_merged_pod5, mod_aligned_bam, mod_chunk_path, can_base,
        mod_symbol=mod_symbol, kmer_context=kmer_context,
        kmer_table_path=kmer_table_path, is_control=False,
        max_chunks_per_read=(2*mod_chunk_range+1)
    )
    prepare_chunks(
        can_merged_pod5, can_aligned_bam, can_chunk_path, can_base,
        kmer_context=kmer_context, kmer_table_path=kmer_table_path,
        is_control=True, max_chunks_per_read=(2*can_chunk_range+1)
    )

    # 7) Config + Train
    train_config = make_train_config(can_chunk_path, mod_chunk_path, chunk_dir, weights=(1,1))
    model_path, validation_log_path = train_model_from_config(
        train_config, model_dir, ml_model_path, batch_size=128, device=0
    )
    print("Xemora [DONE] - Model:", model_path)
    print("Xemora [DONE] - Validation log:", validation_log_path)

if __name__ == "__main__":
    main()

