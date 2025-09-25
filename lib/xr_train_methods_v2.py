########################################################################
########################################################################
"""
xr_train_methods.py 

Title: Unpublished work

Perform single-context Xemora model training. This script is method-based 
and uses Dorado for the initial basecalling.

8/23 - first draft of script done, need to test on a different computer. 
To do: add confusion matrices and loss function plotting scripts as outputs (currently there is a temporary confmat.py call, i think both of these graphs should be in 1 new script) 

By: H. Kawabe, N. Kaplan, J. Sumabat, G. Loia, J. A. Marchand

Updated: 11/28/23
"""
########################################################################
########################################################################

import os
import glob
import sys
from pathlib import Path
from xr_tools  import *
from xr_params import *

########################################################################
print('Xemora [Status] - Initializing Xemora...')

# Initialize directories and input paths from command-line arguments
working_dir = os.path.expanduser(sys.argv[1])
xna_raw_dir = os.path.expanduser(sys.argv[2])
xna_ref_fasta = os.path.expanduser(sys.argv[3])
dna_raw_dir = os.path.expanduser(sys.argv[4])
dna_ref_fasta = os.path.expanduser(sys.argv[5])



# Create necessary directories if they don't exist
working_dir = check_make_dir(working_dir)
ref_dir = check_make_dir(os.path.join(working_dir,'references'))
chunk_dir = check_make_dir(os.path.join(working_dir,'chunks'))
model_dir = check_make_dir(os.path.join(working_dir, 'model'))
mod_dir = check_make_dir(os.path.join(working_dir,'modified_preprocess'))
mod_pod_dir = check_make_dir(os.path.join(mod_dir,'pod5'))
mod_bam_dir = check_make_dir(os.path.join(mod_dir,'bam'))
can_dir = check_make_dir(os.path.join(working_dir, 'canonical_preprocess'))
can_pod_dir = check_make_dir(os.path.join(can_dir, 'pod5'))
can_bam_dir = check_make_dir(os.path.join(can_dir, 'bam'))
########################################################################

# Methods

# Step 0: FASTA to xFASTA conversion
def xfasta_conversion(fasta_file, ref_dir):
    """
    Converts an XNA-containing FASTA file to an xFASTA file format using the 
    xr_fasta2x_rc.py script.

    Parameters: 
    fasta_file (str): File path to the XNA-containing FASTA file.
    ref_dir (str): Directory where the xFASTA file will be saved.

    Returns:
    str: File path to the generated xFASTA file.
    """
    if os.path.isfile(os.path.expanduser(fasta_file)):
        xfasta_filepath = os.path.join(ref_dir, 'x'+os.path.basename(fasta_file))
        cmd = 'python lib/xr_fasta2x_rc.py '+os.path.expanduser(fasta_file)+' '+xfasta_filepath
        os.system(cmd)
        return xfasta_filepath
    else: 
        print('Xemora  [ERROR] - Reference fasta xna file not found. Please check if the file exists or verify the file path.')
        sys.exit()

def validate_read_directory(raw_dir):
    """
    Validates the input directory to ensure it contains only FAST5 or POD5 files.

    Parameters:
    raw_dir (str): Directory containing FAST5 or POD5 files.

    Returns: 
    str: The file type present in the directory (either 'fast5' or 'pod5').
    """
    directory_exists = os.path.isdir(raw_dir)
    homogenous_files = True
    filetype = ""
    if directory_exists:
        directory_files = os.listdir(raw_dir)
        ext_list = []
        for file in directory_files:
            ext = file.split('.')[-1]
            ext_list.append(ext)
        uniques = list(set(ext_list))
        if (len(uniques) != 1):
            homogenous_files = False
            print('Xemora [STATUS] - Passed reads directory not homogeneous. File types found: {}'.format(uniques))
        else:
            filetype = uniques[0]
            
    return filetype
    
# Fast5 to pod5 conversion
def cod5_to_fast5(fast5_input, pod_dir, overwrite_pod):
    """
    Converts a directory of FAST5 files to a single POD5 file using the pod5 package.

    Parameters: 
    fast5_input (str): Directory containing FAST5 raw data.
    pod_dir (str): Directory for outputted single POD5 file.
    overwrite_pod (bool): Parameter to allow overwriting the merged POD5 file.

    Returns:
    str: File path to the merged POD5 file.
    """
    pod5_output = os.path.join(pod_dir, os.path.basename(fast5_input)+'.pod5')
    if overwrite_pod or not os.path.exists(pod5_output):
        cmd = 'pod5 convert fast5 --force-overwrite '+fast5_input+'/*.fast5 -o '+pod5_output
        os.system(cmd)
    else:
        print('Xemora [STATUS] - Skipping FAST5 conversion')
    return pod5_output

# Merge pod5 files
def pod5_merge(pod5_input, pod_dir, overwrite_pod):
    """
    Merges multiple POD5 files into a single POD5 file using the pod5 package.

    Parameters: 
    pod5_input (str): Directory containing POD5 raw data.
    pod_dir (str): Directory for outputted single POD5 file.
    overwrite_pod (bool): Parameter to allow overwriting the merged POD5 file.

    Returns:
    str: File path to the merged POD5 file.
    """
    merged_pod5 = os.path.join(pod_dir, os.path.basename(pod5_input)+'.pod5')
    if overwrite_pod or not os.path.exists(merged_pod5):
        cmd = 'pod5 merge --force-overwrite '+ pod5_input+'/*.pod5 -o ' + merged_pod5
        os.system(cmd)
    else:
        print('Xemora [STATUS]- Skipping POD5 merge')
    return merged_pod5
    
def dorado_basecall(dorado_path, dorado_model, min_qscore, pod_dir, bam_directory, basecall_pod, max_reads, filter_readIDs):
    """
    dorado_basecall takes in various inputs needed to run the Dorado basecaller
    and returns the file path to the generated BAM file 
    
    Parameters: 
    dorado_path - file pathway to Dorado
    dorado_model - canonical model to use for basecalling, refer to xr_params for more info
    xfasta_path - file pathway to an xFASTA file 
    min_qscore - minimum read quality score to allow in BAM file, set in xr_params
    pod_dir - inputted POD5 file or directory to basecall 
    bam_directory - desired output directory for outputted BAM file. 
    basecall_pod - parameter allowing for basecalling to be repeated
    max_reads - maxium number of reads to basecall, default all. In xr_params
    filter_readIDs -  file path to new line delimited readIDs, default: all

    Returns:
    output_bam - BAM file output file pathway
    """
    output_bam = os.path.join(bam_directory, 'bc.bam')
    
    if basecall_pod or not os.path.exists(output_bam):
        print('Xemora [STATUS] Performing basecalling using Dorado')
        #Base arguments
        dorado_args = f'--no-trim --emit-moves  {pod_dir}'

        #Optional arguments
        if filter_readIDs: 
            dorado_args +=f' -l {filter_readIDs}'
        if max_reads:
            dorado_args +=f' -n {max_reads}'
        if min_qscore: 
            dorado_args +=f' --min-qscore {min_qscore}'

        #Dorado command 
        cmd = f'{dorado_path} basecaller {dorado_model} {dorado_args} > {output_bam}'
        os.system(cmd) 
        return output_bam
        
    else:
        print('Xemora [STATUS] - Skipping POD5 basecalling for modified bases.')
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
    
def ensure_index(bam_path):
    bai = bam_path + ".bai"
    csi = bam_path + ".csi"
    if not (os.path.exists(bai) or os.path.exists(csi)):
        print(f"[INFO] Index missing, indexing {bam_path}")
        os.system(f"samtools index {bam_path}")
    else:
        print(f"[INFO] Index present for {bam_path}")


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
    
from typing import Optional

def rc_xfasta_candidate(path: str) -> Optional[str]:



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
        print('Xemora  [STATUS] - Generating chunks for modified basecalling.')
        cmd = 'remora \
          dataset prepare \
          '+pod_file+' \
          '+bam_file+' \
          --output-remora-training-file '+chunk_file+' \
          --focus-reference-positions '+bed_file+' \
          --mod-base '+mod_base+' '+mod_base+' \
          --kmer-context-bases '+kmer_context+' \
          --max-chunks-per-read '+str(2*mod_chunk_range+1)+' \
          --refine-kmer-level-table '+kmer_table_path+' \
          --refine-rough-rescale '+' \
          --motif '+can_base+' 0 \
          --chunk-context '+chunk_context
        os.system(cmd)
        return chunk_file
    else:
        print('Xemora [STATUS] - Skipping modified chunk generation')
        return chunk_file

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
        print('Xemora [STATUS] - Generating chunks for canonical basecalling.') 
        cmd = 'remora \
            dataset prepare \
            '+pod_file+' \
            '+bam_file+' \
            --output-remora-training-file '+chunk_file+' \
            --focus-reference-positions '+bed_file+' \
            --max-chunks-per-read '+str(2*mod_chunk_range+1)+' \
            --mod-base-control \
            --motif '+can_base+' 0 \
            --kmer-context-bases '+kmer_context+' \
            --refine-kmer-level-table '+kmer_table_path+' \
            --refine-rough-rescale '+' \
            --chunk-context '+chunk_context
        os.system(cmd)
    else:
        print('Xemora [STATUS] - Skipping canonical chunk generation')
    
    return chunk_file


def merge_chunks(chunk_dir, mod_chunks, can_chunks, balance_chunks):
    """
    Merges modified and canonical chunks into a single training dataset.

    Parameters:
    chunk_dir (str): Directory containing chunk files.
    mod_chunks (str): File path to the modified chunks.
    can_chunks (str): File path to the canonical chunks.
    balance_chunks (bool): Parameter to balance chunks during merging.

    Returns:
    str: File path to the merged training chunk dataset.
    """
    training_chunks = os.path.join(chunk_dir, 'training_chunks.npz')
    if remerge_chunks:
        print('Xemora [STATUS] - Merging chunks for training.')
        if balance_chunks:
            cmd = 'remora \
              dataset merge \
              --balance \
              --input-dataset '+os.path.join(chunk_dir,'mod_chunks.npz')+' '+chunk_num+'_000 \
              --input-dataset '+os.path.join(chunk_dir,'can_chunks.npz')+' '+chunk_num+'_000 \
              --output-dataset '+os.path.join(chunk_dir,'training_chunks.npz')
            os.system(cmd)
        else:
            cmd = 'remora \
              dataset merge \
              --input-dataset '+mod_chunks+' '+chunk_num+'_000 \
              --input-dataset '+can_chunks+' '+chunk_num+'_000 \
              --output-dataset '+training_chunks
            os.system(cmd)
    else: 
        print('Xemora [STATUS] - Skipping modified and canonical chunk merging')
    return training_chunks

def xemora_training(model_dir, training_chunks):
    """
    Trains the Xemora model using the merged training chunks.

    Parameters:
    model_dir (str): Directory to save the trained model.
    training_chunks (str): File path to the merged training chunk dataset.

    Returns: 
    renamed_model_path - file path to generated single context model 
    val_log_path - file path to validation.log containing training information
    """
    if gen_model:
        print('Xemora [STATUS] - Training model.')
        cmd = 'remora \
          model train \
          '+os.path.join(chunk_dir,'training_chunks.npz')+' \
          --model '+ml_model_path+' \
          --device 0 \
          --output-path '+model_dir+' \
          --overwrite \
          --kmer-context-bases '+kmer_context+' \
          --chunk-context '+chunk_context+' \
          --val-prop '+val_proportion+' \
          --batch-size 100 '# + '

        os.system(cmd)
    else: 
        print('Xemora [STATUS] - Skipping model training')

    model_path =os.path.join(model_dir, 'model_best.pt')
    #renamed_model_path = os.path.join(model_dir, mod_base+can_base+'_model.pt')
    #os.rename(model_path, renamed_model_path)
    val_log_path = os.path.join(model_dir, 'validation.log')

    return model_path, val_log_path


def run_remora_ref_region_plot(
    can_pod5, can_bam,
    mod_pod5, mod_bam,
    ref_bed,
    highlight_bed,
    levels_table,
    out_dir,
    prefix="ref_region",
    log_name="ref_region.log"
):
    """
    Generate Remora ref_region plots (Remora 2.1 style).

    Parameters
    ----------
    can_pod5 : str
        Path to canonical pod5 file
    can_bam : str
        Path to canonical aligned bam file
    mod_pod5 : str
        Path to modified pod5 file
    mod_bam : str
        Path to modified aligned bam file
    ref_bed : str
        BED file with ref regions to plot
    highlight_bed : str
        BED file with highlight regions (optional, can be same as ref_bed)
    levels_table : str
        Path to k-mer levels table
    out_dir : str
        Directory to save plots
    prefix : str
        Prefix for output files
    log_name : str
        Name of log file
    """
    os.makedirs(out_dir, exist_ok=True)
    cwd = os.getcwd()
    try:
        # Resolve absolute paths
        ref_bed = os.path.abspath(os.path.expanduser(ref_bed))
        highlight_bed = os.path.abspath(os.path.expanduser(highlight_bed))
        levels_table = os.path.abspath(os.path.expanduser(levels_table))

        if not os.path.exists(levels_table):
            raise FileNotFoundError(f"[ERROR] k-mer levels table not found: {levels_table}")

        os.chdir(out_dir)

        # Build the Remora command
        cmd = (
            f"remora analyze plot ref_region "
            f"--pod5-and-bam {can_pod5} {can_bam} "
            f"--pod5-and-bam {mod_pod5} {mod_bam} "
            f"--ref-regions {ref_bed} "
            f"--highlight-ranges {highlight_bed} "
            f"--refine-kmer-level-table {levels_table} "
            f"--refine-rough-rescale "
            f"--log-filename {log_name}"
        )
        print(f"[DEBUG] Running Remora ref_region plot:\n{cmd}")
        rc = os.system(cmd)
        if rc != 0:
            raise RuntimeError(f"Remora ref_region plotting failed with code {rc}")
        print(f"Xemora [STATUS] - Remora ref_region plots written to {out_dir}")
    finally:
        os.chdir(cwd)



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

    ensure_index(mod_aligned_bam)
    ensure_index(can_aligned_bam)

        # Step: Generate Remora plots
    remora_plot_dir = check_make_dir(os.path.join(working_dir, "remora_plots"))
    run_remora_ref_region_plot(
        can_merged_pod5, can_aligned_bam,
        mod_merged_pod5, mod_aligned_bam,
        ref_bed=ref_regions_bed,             # or mod_bed_file if you want mod-specific
        highlight_bed=mod_bed_file,       # whichever highlight you want
        levels_table=kmer_table_path,
        out_dir=remora_plot_dir
    )

    # Step 4: Generate Remora chunks
    mod_chunk_path = generate_mod_chunks(mod_merged_pod5, mod_aligned_bam, chunk_dir, mod_bed_file, mod_base, kmer_context, kmer_table_path, regenerate_chunks)
    can_chunk_path = generate_can_chunks(can_merged_pod5, can_aligned_bam, chunk_dir, can_bed_file, can_base, kmer_context, kmer_table_path, regenerate_chunks)


    # Step 5: Merge chunk files
    training_chunk_path = merge_chunks(chunk_dir, mod_chunk_path, can_chunk_path, balance_chunks)

    # Step 6: Train Xemora model
    model_path, validation_log_path = xemora_training(model_dir, training_chunk_path)




if __name__ == "__main__":
    main()

