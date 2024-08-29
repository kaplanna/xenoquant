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
    
def dorado_basecall(dorado_path, dorado_model, xfasta_path, min_qscore, pod_dir, bam_directory, basecall_pod):
    """
    Runs the Dorado basecaller to generate a BAM file from POD5 files.

    Parameters: 
    dorado_path (str): File path to the Dorado executable.
    dorado_model (str): Canonical model to use for basecalling.
    xfasta_path (str): File path to an xFASTA file.
    min_qscore (int): Minimum read quality score to include in the BAM file.
    pod_dir (str): Directory containing the POD5 file(s) to be basecalled.
    bam_directory (str): Desired output directory for the BAM file.
    basecall_pod (bool): Parameter to allow basecalling to be repeated.

    Returns:
    str: File path to the output BAM file.
    """
    output_bam = os.path.join(bam_directory, 'bc.bam')
    
    if basecall_pod or not os.path.exists(output_bam):
        cmd = '{} basecaller {} --reference {} --no-trim --emit-moves --min-qscore {} {} > {}'.format(dorado_path, dorado_model, xfasta_path, min_qscore, pod_dir, output_bam)
        os.system(cmd) 
        return output_bam
    else:
        print('Xemora [STATUS] - Skipping POD5 basecalling for modified bases.')
        return output_bam

# Step 4: Bed file generation 
def bed_gen(input_fasta, xna_base, sub_base, chunk_range, chunk_shift): 
    """
    Generates a BED file specifying the regions to focus on during chunk generation.

    Parameters:
    input_fasta (str): xFASTA file path.
    xna_base (str): XNA of interest, defined in xr_params.
    sub_base (str): Base to substitute with the XNA base.
    chunk_range (int): Integer representing the +/- range around focus to analyze.
    chunk_shift (int): Integer representing how many bases to shift the focus region.

    Returns:
    str: File path to the generated BED file.
    """
    output_bed = os.path.join(os.path.dirname(input_fasta), xna_base+'.bed')
    if os.stat(input_fasta).st_size == 0: 
        print('Xemora  [ERROR] - Empty xfasta file generated. Check that XNA bases were present in sequence of input fasta file.')
        sys.exit()
    else:
        fr = open(output_bed,"w")
        with open(os.path.expanduser(input_fasta), "r") as fo:
            for line in fo: 
                if 'GAP' not in line.upper(): 
                    if line[0]=='>':
                        header = line[1:].replace('\n','')
                        x_pos_base = fetch_xna_pos(header)
                        x_pos_to_rc =[]

                        for x in x_pos_base: 
                            x_base = x[0]
                            x_pos = int(''.join(filter(str.isdigit, x[1])))

                            if x_base == xna_base: 
                                strand = '+'
                            elif x_base == xna_base_rc(xna_base,xna_base_pairs): 
                                strand ='-'
                            bed_line = header+'\t'+str(x_pos-chunk_range+chunk_shift)+'\t'+str(int(x_pos)+chunk_range+1+chunk_shift)+'\t'+sub_base+'\t0\t'+strand+'\n'
                            fr.write(bed_line)

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
    chunk_file = os.path.join(chunk_dir, 'modified_chunks.npz')
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
          --motif N 0 \
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
    chunk_file = os.path.join(chunk_dir, 'canonical_chunks.npz')
    if regenerate_chunks:
        print('Xemora [STATUS] - Generating chunks for canonical basecalling.') 
        cmd = 'remora \
            dataset prepare \
            '+pod_file+' \
            '+bam_file+' \
            --output-remora-training-file '+chunk_file+' \
            --focus_reference-positions '+bed_file+' \
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
          '+training_chunks+' \
          --model '+ml_model_path+' \
          --device 0 \
          --output-path '+model_dir+' \
          --overwrite \
          --kmer-context-bases '+kmer_context+' \
          --chunk-context '+chunk_context+' \
          --val-prop '+val_proportion+' \
          --batch-size 100 '# + '\
    else: 
        print('Xemora [STATUS] - Skipping model training')

    model_path =os.path.join(model_dir, 'model_best.pt')
    renamed_model_path = os.path.join(model_dir, mod_base+can_base+'_model.pt')
    os.rename(model_path, renamed_model_path)
    val_log_path = os.path.join(model_dir, 'validation.log')

    return renamed_model_path, val_log_path

def gen_confmat(val_log_path, output_dir):
    """
    Generates confusion matrices from the validation log.

    Parameters:
    val_log_path (str): File path to the validation log.
    output_dir (str): Directory to save the confusion matrices.
    """
    cmd = 'python xr_confmat.py '+val_log_path+' '+output_dir
    os.system(cmd)

def main():
    '''
    Runs the entire pipeline to generate a reference-anchored (single-context) Xemora model.
    '''
    # Step 0: xFASTA conversion 
    mod_xfasta = xfasta_conversion(xna_ref_fasta, ref_dir) # Modified dataset xFASTA file 
    can_xfasta = xfasta_conversion(dna_ref_fasta, ref_dir) # Canonical dataset xFASTA file

    # Step 1: Merge POD5 files 
    
    # Verify only one file type in raw data directories 
    mod_filetype = validate_read_directory(xna_raw_dir)
    can_filetype = validate_read_directory(dna_raw_dir)

    filetypes = [mod_filetype, can_filetype]

    # Perform either FAST5 conversion or POD5 merging

    # For modified dataset
    if mod_filetype == 'fast5':
        mod_merged_pod5 = cod5_to_fast5(xna_raw_dir, mod_pod_dir, overwrite_pod)
    elif mod_filetype == 'pod5':
        mod_merged_pod5 = pod5_merge(xna_raw_dir, mod_pod_dir, overwrite_pod)
    else:
        print('Xemora [ERROR] - filetype in XNA raw data directory is not POD5 or FAST5, please check if your directory is correct') 
        sys.exit()

    # For canonical dataset
    if can_filetype == 'fast5':
        can_merged_pod5 = cod5_to_fast5(dna_raw_dir, can_pod_dir, overwrite_pod)
    elif can_filetype == 'pod5':
        can_merged_pod5 = pod5_merge(dna_raw_dir, can_pod_dir, overwrite_pod)
    else:
        print('Xemora [ERROR] - filetype in DNA raw data directory is not POD5 or FAST5, please check if your directory is correct') 
        sys.exit()

    # Step 2: Perform basecalling using Dorado
    mod_bc_bam = dorado_basecall(dorado_path, dorado_model, mod_xfasta, min_qscore, mod_merged_pod5, mod_bam_dir, basecall_pod) # Modified dataset
    can_bc_bam = dorado_basecall(dorado_path, dorado_model, can_xfasta, min_qscore, can_merged_pod5, can_bam_dir, basecall_pod)

    # Step 3: Generate BED file for region to analyze 
    mod_bed_file = bed_gen(mod_xfasta, mod_base, mod_base, mod_chunk_range, mod_chunk_shift)
    can_bed_file = bed_gen(can_xfasta, mod_base, can_base, can_chunk_range, can_chunk_shift) # Can shift and range variables don't exist yet, ask Jorge if useful for single-context training 

    # Step 4: Generate Remora chunks
    mod_chunk_path = generate_mod_chunks(mod_merged_pod5, mod_bc_bam, chunk_dir, mod_bed_file, mod_base, kmer_context, kmer_table_path, regenerate_chunks)
    can_chunk_path = generate_can_chunks(can_merged_pod5, can_bc_bam, chunk_dir, can_bed_file, can_base, kmer_context, kmer_table_path, regenerate_chunks)

    # Step 5: Merge chunk files
    training_chunk_path = merge_chunks(chunk_dir, mod_chunk_path, can_chunk_path, balance_chunks)

    # Step 6: Train Xemora model
    model_path, validation_log_path = xemora_training(model_dir, training_chunk_path)

    # Step 7 (Optional): Generate Confusion Matrices from validation.log
    if generate_confusion_matrices:
        confmat_dir = check_make_dir(os.path.join(working_dir, 'confusion_matrices'))
        gen_confmat(validation_log_path, confmat_dir)

