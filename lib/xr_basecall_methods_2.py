########################################################################
########################################################################
"""
xr_basecall_pairwise.py 

Title: Unpublished work

Perform 'standard' single context basecalling for multiple model inputs. Also handles 
pairwise analysis if desired. 


Moved logs from old xr_basecall_methods
[Change log] 
9/4 - added chck_fasta, copy_fasta, and full_read_bed functions to either perform full sequence inference or hard swap
9/4 - during chunk generation, removed dynamic range analysis and set to 1000 (since the range is not constant for different reads) 
9/4 - doc strings added to new functions for decoupled consensus analysis 
9/4 - general print statements added/edited for more clarity on what is going on
9/4 - removed focus base list and dataframe printing from consensus analysis, not useful
9/4 - moved min_cluster_size to xr_params
9/5 - min_cluster_size renamed to min_coverage 
9/5 - mod_chunk_range is now always used. For full read in bed file , set mod_chunk_range = -1 (only for non-xfasta references)
9/5 - Added minimap filters: min_mapq_score and min_as_score. Both are used to filter the bam output. Can be set in xr_params
9/5 - changed variable name of 'xna_ref_fasta' to just ref_fasta to better match new added functionality
9/17 - reverse strand handling fixed for bed file chunk shift, confirmed with indexing math and remora focus_base outputs 
9/17 - Shannon Entropy and QS extraction fixed to obtain correct indexes, confirmed by mode of extracted bases compared to reference file
9/17 - locking this copy as working, begining development for pairwise analysis in xr_basecall_pairwise
9/25 - pairwise logic fixedd 
9/26 edited xemora.py to handle multiple model file inputs

[To Do] 
- Change directory structure to mirror basecall methods: chunk, preprocessing, remora results (all models), consensus results (all models).
- Add flag to do pair wise that only works if 2 models are inputted 
- bed file generation generalization 

-- if no alignments in bam file post mm2, have it exit script
[Logs]
9/17 - file generated

By: H. Kawabe, N. Kaplan, J. Sumabat, J. A. Marchand

Updated: 11/28/23
"""
########################################################################
########################################################################

import os
import glob
import sys
import json
import pysam
import sys
import tempfile
from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path
from xr_tools  import *
from xr_params import *
from remora.data_chunks import RemoraRead
from remora.model_util import load_model
from remora.inference import call_read_mods
from collections import defaultdict


########################################################################

#Basecall paths
working_dir = os.path.expanduser(sys.argv[1])
raw_dir = os.path.expanduser(sys.argv[2])
ref_fasta = os.path.expanduser(sys.argv[3])
model_files = sys.argv[4:-1]  # Take in all model files from position 4 until params
params = json.loads(sys.argv[-1]) if len(sys.argv) > 4 else {}


############################################################

#########################################################################

#Generate directories
print('Xemora [Status] - Initializing Xemora basecalling.')

working_dir = check_make_dir(working_dir)
ref_dir = check_make_dir(os.path.join(working_dir,'references'))
preprocess_dir = check_make_dir(os.path.join(working_dir,'preprocess'))
pod_dir = check_make_dir(os.path.join(preprocess_dir,'pod5'))
mod_bam_dir = check_make_dir(os.path.join(preprocess_dir, 'bam'))
chunk_dir = check_make_dir(os.path.join(working_dir, 'chunks')) 
remora_output_dir = check_make_dir(os.path.join(working_dir, 'remora_outputs'))
analysis_dir = check_make_dir(os.path.join(working_dir, 'results'))

#########################################################################
#Methods 
def set_params(params='{}'):
    """
    set_params lets you overwrite parameters from xr_params by inputting a dict
    or a JSON string containing the parameters you wish to overwrite.
    """
    if isinstance(params, str):
        try:
            params_dict = json.loads(params)
        except json.JSONDecodeError as e:
            print(f"Xemora [ERROR] - Invalid JSON string for parameters: {params}")
            sys.exit(1)
    elif isinstance(params, dict):
        params_dict = params
    else:
        print(f"Xemora [ERROR] - Unsupported type for parameters: {type(params)}")
        sys.exit(1)

    overwritten_params = []
    for key, value in params_dict.items():
        globals()[key] = value
        overwritten_params.append(f"{key}={value}")
    
    if overwritten_params:
        print("Xemora [Status] - Overwritten parameters: " + ", ".join(overwritten_params))
    else:
        print("Xemora [Status] - Default parameter sets from xr_params used.")


#check fasta file if XNAs are present in the sequences 
def check_fasta(fasta_file, xna_base_pairs):
    """
    check_fasta looks at the inputted FASTA file and checks for the presence of 
    XNAs based on the list provided in xna_base_pairs. Each pair represents
    a non-standard base, and the function returns True if any XNA is found in the sequences.

    Parameters:
    fasta_file (str): Path to the input FASTA file.
    xna_base_pairs (list of str): List of XNA base pairs, where each pair represents a combination of non-standard nucleobases. from xr_params

    Returns:
    bool: True if an XNA base is found in any of the sequences, False otherwise.
    """
    
    print('Xemora [STATUS] - Checking inputted fasta file for XNAs')
    # Convert the xna_base_pairs to a set of individual XNA characters
    xna_bases = set(''.join(xna_base_pairs))
    
    # Parse the FASTA file
    with open(fasta_file, "r") as fh:
        for record in SeqIO.parse(fh, "fasta"):
            sequence = set(str(record.seq).upper())  # Convert sequence to a set of unique bases
            
            # Check if any XNA base is present in the sequence
            if xna_bases.intersection(sequence):
                print('Xemora [STATUS] - XNA found in input fasta, performing hardswap')
                return True  # Return True if there is any intersection (XNA found)
    
    # If no XNA is found after checking all sequences
    print('Xemora [STATUS] - No XNA found in input fasta, performing all base validation')
    return False

#Step 0: FASTA to xFASTA conversion
def xfasta_conversion(fasta_file, ref_dir):
    """
    xfasta_conversion takes an XNA containing FASTA file and generates an 
    xFASTA formatted fil by calling the xr_fasta2x_rc.py script
    
    Parameters: 
    fasta_file - file pathway to XNA containing fasta file 
    ref_dir - desired reference/output directory for the xFASTA file 
    
    Returns:
    xfasta_filepath - filepath to generated xFASTA file.
    """
    if os.path.isfile(os.path.expanduser(fasta_file)):
    
        xfasta_filepath = os.path.join(ref_dir, 'x'+os.path.basename(fasta_file))
        cmd = 'python lib/xr_fasta2x_rc.py '+os.path.expanduser(fasta_file)+' '+xfasta_filepath
        os.system(cmd)
        return xfasta_filepath
    else: 
        print('Xemora  [ERROR] - Reference fasta xna file not file. Please check file exist or file path.')
        sys.exit()

def copy_fasta(ref_fasta, ref_dir):
    """
    Copies the inputted FASTA file to the ref_dir, generates a reverse complement
    version, and returns the file path along with a dictionary of sequence names and lengths.
    
    Parameters:
    ref_fasta (str): Path to the input FASTA file.
    ref_dir (str): Directory where the copied and reverse complement files will be saved.
    
    Returns:
    copied_fasta_path (str): Path to the copied FASTA file.
    ref_lengths (dict): Dictionary where keys are sequence names and values are sequence lengths.
    """
    
    # Ensure the directory exists
    if not os.path.exists(ref_dir):
        os.makedirs(ref_dir)
    
    copied_fasta_path = os.path.join(ref_dir, os.path.basename(ref_fasta))
    rc_fasta_path = os.path.join(ref_dir, 'rc_' + os.path.basename(ref_fasta))
    
    ref_lengths = {}
    
    with open(copied_fasta_path, 'w') as copy_fh, open(rc_fasta_path, 'w') as rc_fh:
        # Parse the input FASTA file and process sequences
        for record in SeqIO.parse(ref_fasta, "fasta"):
            # Write the original sequence
            SeqIO.write(record, copy_fh, "fasta")
            
            # Calculate and store the length of the sequence
            ref_lengths[record.id] = len(record.seq)
            
            # Create reverse complement sequence
            rc_seq = record.seq.reverse_complement()
            rc_record = record
            rc_record.seq = rc_seq
            rc_record.id = record.id + "_rc"  # Append '_rc' to the ID for reverse complement
            
            # Write the reverse complement sequence
            SeqIO.write(rc_record, rc_fh, "fasta")
    
    return copied_fasta_path, ref_lengths

def validate_read_directory(raw_dir):
    """
    validate_read_directory takes in the FAST5 or POD5 directory input and verifies
    only FAST5 and POD5 data is present. This method also outputs the file type 
    that is present in the directory
    
    Parameters:
    raw_dir - directory containing FAST5 or POD5 files 
    
    Returns: 
    filetype - filetype present 
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
            print('Xemora [STATUS] - Passed reads directory not homogenous. Filetypes found: {}'.format(uniques))
        else:
            filetype = uniques[0]
            
    return filetype
    
#Fast5 to pod5 conversion
def cod5_to_fast5(fast5_input, pod_dir, overwrite_pod):
    """
    cod5_to_fast5 takes in a FAST5 directory and a desired output directory and 
    runs the pod5 convert fast5 command from the pod5 package 
    
    Parameters: 
    fast5_input - directory containing FAST5 raw data 
    pod_dir - directory for outputted single POD5 file 
    overwrite_pod - parameter from xr_params to allow the merged POD5 file to be overwritten 
    
    Returns:
    pod5_output - file pathway to merged POD5 file 
    """
    pod5_output = os.path.join(pod_dir, os.path.basename(fast5_input)+'.pod5')
    if overwrite_pod or not os.path.exists(pod5_output):
        cmd = 'pod5 convert fast5 --force-overwrite '+fast5_input+'/*.fast5 -o '+pod5_output
        os.system(cmd)
    else:
        print('Xemora [STATUS] - Skipping FAST5 conversion')
    return pod5_output

#Merge pod5 files
def pod5_merge(pod5_input, pod_dir, overwrite_pod):
    """
    pod5_merge takes in a FAST5 directory and a desired output directory and runs 
    the pod5 merge command from the pod5 package 
    
    Parameters: 
    pod5_input - directory containing POD5 raw data 
    pod_dir - directory for outputted single POD5 file 
    overwrite_pod - parameter from xr_params to allow the merged POD5 file to be overwritten
    """
    merged_pod5 = os.path.join(pod_dir, os.path.basename(pod5_input)+'.pod5')
    if overwrite_pod or not os.path.exists(merged_pod5):
        cmd = 'pod5 merge --force-overwrite '+ pod5_input+'/*.pod5 -o ' + merged_pod5
        os.system(cmd)
    else:
        print('Xemora [STATUS] - Skipping POD5 merge')
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
        # Add Dorado summary step here after basecalling
        summary_cmd = 'dorado summary -v {} > {}'.format(output_bam, os.path.join(bam_directory, 'sequencing_summary.txt'))
        os.system(summary_cmd)
        print('Xemora [STATUS] - Dorado summary saved to sequencing_summary.txt.')
        
        return output_bam
        
    else:
        print('Xemora [STATUS] - Skipping POD5 basecalling for modified bases.')
        return output_bam




def filter_bam_by_qs(bam_file_path, output_bam_directory):
    """
    filter_bam_by_qs takes in a BAM file pathway, a desired output directory 
    and a quality score threshold and filters out reads with a lower mean quality
    score 

    Parameters:
    bam_file_path - Inputted BAM file path as a string 
    output_bam_directory - Desired output directory for filtered BAM file 
    min_qscore_analysis - Post-basecalling filter of reads by q-score

    Returns:
    qs_filered_bam - output BAM file pathway
    """
    qs_filtered_bam = os.path.join(output_bam_directory, 'qs_filtered.bam')
    # Open input BAM file
    with pysam.AlignmentFile(bam_file_path, "rb", check_sq=False) as bam_in:
        # Create output BAM file
        with pysam.AlignmentFile(qs_filtered_bam, "wb", header=bam_in.header) as bam_out:
            for read in bam_in:
                # Extract 'qs:f' tag (mean basecall qscore)
                if read.has_tag("qs:f"):
                    qscore = read.get_tag("qs:f")
                    # Write to output if qscore is greater than the threshold
                    if qscore >= min_qscore_analysis:
                        bam_out.write(read)
    return qs_filtered_bam

def minimap2_aligner(input_bam, fasta_path, bam_directory):
    """
    minimap2_aligner takes in a bam file from Dorado basecaller and outputs a new 
    aligned bam file using an xFASTA formatted file 

    Parameters:
    input_bam - file path to bam file to align
    xfasta_path - file pathway to FASTA file to align basecalls to 
    bam_directory - output file pathway for aligned BAM file 

    Returns: 
    aligned_bam - minimap2 aligned bam file 
    """
    print('Xemora [STATUS] - Aligning basecalls using Minimap2')
    aligned_bam = os.path.join(bam_directory, 'aligned.bam')
    cmd = f'samtools fastq -T "*" {input_bam} | minimap2 -y -ax map-ont --score-N 0 --secondary no --sam-hit-only --MD {fasta_path} - | samtools view -q {min_mapq_score} -F0x800 -bho {aligned_bam}'
    #cmd = f'samtools fastq -T "*" {input_bam} | minimap2 -y -ax map-ont --score-N 0 --secondary no --sam-hit-only --MD {xfasta_path} - | samtools view -F0x800 -bho {aligned_bam}'
    #cmd = 'samtools fastq -T "*" '+input_bam+ ' | minimap2 -y -ax map-ont --score-N 0  --MD '+xfasta_path+ ' - | samtools view -bho ' +aligned_bam #checking for unaligned reads generated 
    os.system(cmd)

    return aligned_bam

def filter_by_alignment_score(bam_file):
    """
    Filter reads based on the alignment score (AS) field and rewrite the original BAM file.

    Parameters:
        bam_file (str): Path to the input BAM file.
        min_alignment_score (int): Minimum alignment score to keep a read.
        
    Raises:
        FileNotFoundError: If the BAM file does not exist.
        ValueError: If the BAM file cannot be parsed.
    """
    # Suppress warnings about missing BAM index
    pysam.set_verbosity(0)
    
    try:
        # Open the BAM file
        bam = pysam.AlignmentFile(bam_file, "rb")
    except FileNotFoundError:
        raise FileNotFoundError(f"BAM file {bam_file} does not exist.")
    except ValueError:
        raise ValueError(f"Unable to parse the BAM file {bam_file}.")

    # Create a temporary file for the filtered BAM
    temp_bam_fd, temp_bam_path = tempfile.mkstemp(suffix=".bam")
    os.close(temp_bam_fd)

    try:
        # Open the temporary BAM file for writing
        with pysam.AlignmentFile(temp_bam_path, "wb", header=bam.header) as out_bam:
            # Iterate through each read in the BAM file
            for read in bam:
                # Check if the read is a primary alignment
                if not read.is_secondary and not read.is_supplementary:
                    # Get the alignment score (AS) field
                    alignment_score = read.get_tag('AS') if read.has_tag('AS') else None
                    
                    # Check if the alignment score meets the minimum threshold
                    if alignment_score is not None and alignment_score >= min_as_score:
                        out_bam.write(read)

        # Close the BAM file
        bam.close()

        # Replace the original BAM file with the filtered BAM file
        os.replace(temp_bam_path, bam_file)

    finally:
        # Clean up the temporary file if it still exists
        if os.path.exists(temp_bam_path):
            os.remove(temp_bam_path)

    # Reset verbosity
    pysam.set_verbosity(1)

def get_primary_alignments(bam_file):
    """
    Extracts primary alignments from a BAM file, including reference lengths.

    This function opens a BAM file and iterates through each read, 
    identifying the primary alignments (i.e., alignments that are not 
    secondary or supplementary). It maps the read IDs to a list containing 
    its primary reference sequence (RNAME), CIGAR string, read sequence, 
    quality scores, alignment position (POS), and the reference length.

    Parameters:
        bam_file (str): Path to the BAM file.

    Returns:
        dict: A dictionary where the keys are read IDs (query names) and the 
        values are lists containing the reference sequence name (RNAME), 
        CIGAR string, read sequence, quality scores, the alignment position (POS), 
        and the reference length (REF_LEN).
        
    Raises:
        FileNotFoundError: If the BAM file does not exist.
        ValueError: If the BAM file cannot be parsed.
    """
    # Suppress warnings about missing BAM index
    pysam.set_verbosity(0)
    
    try:
        # Open the BAM file
        bam = pysam.AlignmentFile(bam_file, "rb")
    except FileNotFoundError:
        raise FileNotFoundError(f"BAM file {bam_file} does not exist.")
    except ValueError:
        raise ValueError(f"Unable to parse the BAM file {bam_file}.")

    # Dictionary to store read_id -> [reference sequence, CIGAR string, read sequence, QUAL, POS, REF_LEN]
    read_to_info = {}

    # Iterate through each read in the BAM file
    for read in bam:
        # Check if the read is a primary alignment
        if not read.is_secondary and not read.is_supplementary:
            read_id = read.query_name
            ref_name = bam.get_reference_name(read.reference_id)  # RNAME
            cigar_string = read.cigarstring  # CIGAR
            basecalled_sequence = read.query_sequence  # SEQ
            q_score = read.query_qualities  # QUAL
            ref_start_pos = read.reference_start  # POS (0-based)
            flag = read.flag

            # Get the reference length for the current reference
            ref_length = bam.lengths[read.reference_id]  # REF_LEN

            # Store the read_id and associated information, including the reference length
            read_to_info[read_id] = {
                'reference_sequence': ref_name,
                'cigar_string': cigar_string,
                'basecalled_sequence': basecalled_sequence,
                'q_score': q_score,
                'ref_start_pos': ref_start_pos,
                'flag': flag,
                'reference_length': ref_length  # Add reference length
            }

    # Close the BAM file
    bam.close()

    # Reset verbosity
    pysam.set_verbosity(1)

    return read_to_info

def load_xna_model(model_path):
    """
    load_xna_model uses Remora's API to select the XNA and canonical base to perform 
    analysis on 

    Parameters:
    model_path - file path way to model path way 
    
    Returns 
    mod_base - XNA of interest
    can_base - standard base of interest
    """
    #Load model and use model
    print('Xemora [STATUS] - Loading Xemora model for basecalling...')
    model, model_metadata = load_model(model_path)
    mod_base = model_metadata['mod_bases']
    can_base = model_metadata['motif_0']
    print(f'Xemora [STATUS] - Loaded {model_path} trained on: XNA: {mod_base} / DNA: {can_base}')

    return mod_base, can_base 

def bed_gen(input_fasta, xna_base, sub_base, chunk_range, chunk_shift, model_number): 
    """
    bed_gen takes in an xFASTA file, xna base, a base to substitute with the 
    xna base, region for analysis around focus, and shift from focus to generate 
    BED files with the regions to perform Remora chunk generation on both forward 
    and reverse strands.
    
    Parameters:
    input_fasta - xFASTA file path 
    xna_base - XNA of interest
    sub_base - base to substitute
    chunk_range - integer +/- around focus to analyze 
    chunk_shift - integer with how many bases to move the focus region
    
    Returns:
    bed_paths - list containing 2 file paths, containing forward and reverse bed paths 
    """

    # Define output paths for forward and reverse BED files
    fwd_bed = os.path.join(os.path.dirname(input_fasta), f'{xna_base}_model_{model_number}_fwd.bed')
    rev_bed = os.path.join(os.path.dirname(input_fasta), f'{xna_base}_model_{model_number}_rev.bed')

    bed_paths = [fwd_bed, rev_bed]
    if os.stat(input_fasta).st_size == 0: 
        print('Xemora [ERROR] - Empty xFASTA file generated. Check that XNA bases were present in the sequence of input fasta file.')
        sys.exit()

    if chunk_range < 0: 
        print('Xemora [ERROR] - Current chunk_range < 0. Set chunk_range >= 0 for proper chunk extraction.')
        sys.exit()
    
    with open(fwd_bed, "w") as fwd_fr, open(rev_bed, "w") as rev_fr:
        with open(os.path.expanduser(input_fasta), "r") as fo:
            for line in fo: 
                if 'GAP' not in line.upper(): 
                    if line[0] == '>':
                        header = line[1:].strip()  # Sequence header
                        x_pos_base = fetch_xna_pos(header)  # Fetch XNA positions

                        for x in x_pos_base: 
                            x_base = x[0]  # Base letter
                            x_pos = int(''.join(filter(str.isdigit, x[1])))  # Base position

                            # Write the forward strand BED line
                            fwd_bed_line = f'{header}\t{x_pos - chunk_range + chunk_shift}\t{x_pos + chunk_range + 1 + chunk_shift}\t{sub_base}\t0\t+\n'
                            fwd_fr.write(fwd_bed_line)

                            # Write the reverse strand BED line
                            rev_bed_line = f'{header}\t{x_pos - chunk_range - chunk_shift}\t{x_pos + chunk_range + 1 - chunk_shift}\t{sub_base}\t0\t-\n'
                            rev_fr.write(rev_bed_line)

    return bed_paths

def full_read_bed(ref_lengths, xna_base, sub_base, ref_dir, mod_chunk_range=0, mod_chunk_shift=0):
    """
    Generates two BED files (one for forward strand, one for reverse strand) 
    with one line for each reference sequence. This function assumes all sequences
    are on both strands.

    Parameters:
    ref_lengths (dict): Dictionary with reference names and their lengths.
    xna_base (str): The XNA base of interest.
    sub_base (str): The base to substitute the XNA with.
    ref_dir (str): The directory where the BED files will be saved.
    mod_chunk_range (int, optional): Range around the center for chunking. Defaults to 0.
    mod_chunk_shift (int, optional): Shift applied to the chunk. Defaults to 0.

    Returns:
    bed_paths (list): List containing the paths to the generated forward and reverse strand BED files.
    """

    fwd_bed_path = os.path.join(ref_dir, f'{xna_base}_model_{model_number}_fwd.bed')
    rev_bed_path = os.path.join(ref_dir, f'{xna_base}_model_{model_number}_rev.bed')

    with open(fwd_bed_path, 'w') as fwd_bed_fh, open(rev_bed_path, 'w') as rev_bed_fh:
        for ref_name, length in ref_lengths.items():
            read_center = round(length / 2)

            if mod_chunk_range >= 0:
                # Forward strand BED line
                fwd_bed_line = f"{ref_name}\t{read_center - mod_chunk_range + mod_chunk_shift}\t{read_center + mod_chunk_range + mod_chunk_shift}\t{xna_base}\t{sub_base}\t+\n"
                fwd_bed_fh.write(fwd_bed_line)

                # Reverse strand BED line
                rev_bed_line = f"{ref_name}\t{read_center - mod_chunk_range - mod_chunk_shift}\t{read_center + mod_chunk_range - mod_chunk_shift}\t{xna_base}\t{sub_base}\t-\n"
                rev_bed_fh.write(rev_bed_line)
            else:
                # Write full sequence length for both strands
                fwd_bed_line = f"{ref_name}\t0\t{length}\t{xna_base}\t{sub_base}\t+\n"
                fwd_bed_fh.write(fwd_bed_line)

                rev_bed_line = f"{ref_name}\t0\t{length}\t{xna_base}\t{sub_base}\t-\n"
                rev_bed_fh.write(rev_bed_line)

    # Store the paths to the two BED files in a list
    bed_paths = [fwd_bed_path, rev_bed_path]

    return bed_paths

def generate_chunks(pod_file, bam_file, chunk_file, bed_file, mod_base, kmer_context, kmer_table_path, regenerate_chunks):
    """
    generate_chunks prepares and runs a Remora command for chunk generation.
    
    Parameters:
    pod_file - merged POD5 file pathway
    bam_file - aligned BAM file basecalled from POD5 file 
    bed_file - BED file with region/strand to analyze 
    mod_base - modified base to perform inference on, in xr_params
    kmer_context - in xr_patams 
    kmer_table_path - in xr_params 
    regenerate_chunks - parameter to overwrite generated chunk files 
    
    Returns 
    chunk_file - file pathway for Remora chunks 
    
    To do
    Either a parameter or use xr_params XNA:DNA pairs to perform inference
    Come up with something for N basecalling
    """
    #Step 5: Generate Chunks. 
    if regenerate_chunks == True:
        print('Xemora [STATUS] - Generating chunks for modified basecalling.')
        #no motif
        cmd = 'remora \
          dataset prepare \
          '+pod_file+' \
          '+bam_file+' \
          --output-remora-training-file '+chunk_file+' \
          --focus-reference-positions '+bed_file+' \
          --mod-base '+mod_base+' '+mod_base+' \
          --kmer-context-bases '+kmer_context+' \
          --max-chunks-per-read 1000 \
          --refine-kmer-level-table '+kmer_table_path+' \
          --refine-rough-rescale '+' \
           --motif '+can_base+' 0 \
          --chunk-context '+chunk_context
        os.system(cmd)
        return chunk_file
    else:
        print('Xemora [STATUS] - Skipping chunk generation')
        return chunk_file

def xemora_basecall(chunk_file, model_file, per_read_output, summary_output):
    """
    xemora_basecall takes in the generate chunk file and a Xemora model and performs 
    inference on the chunk file to determine is XNAs are present 
    
    Parameters: 
    remora_output_dir - output directory inference files
    chunk_file - file pathway to a Remora chunk file 
    model_file - file pathway to a Xemora model 
    
    Returns: 
    """
    print('Xemora [STATUS] - Performing basecalling.')
    try:
        cmd = f'remora validate from_remora_dataset {chunk_file} --model {os.path.expanduser(model_file)} --full-results-filename {per_read_output} --out-file {summary_output}'
        os.system(cmd)

        print('Xemora [STATUS] - Basecalling done. Saving read-level results to'+per_read_output)
        print('Xemora [STATUS] - Basecalling done. Saving Remora summary to'+summary_output)
        #print('Xemora [STATUS] - Exiting')
        return per_read_output, summary_output
        
    except: 
        print('Xemora [ERROR] - Failed to initialize basecalling model Check logs.')
        sys.exit()

def add_per_read_mapping(primary_alignments, per_read_modifications):
    """
    add_per_read_mapping adds reference sequences to the per_read_modifications 
    file outputted from Remora. Used to perform decoupled analysis downstream
    """
    # Read the CSV file into a DataFrame
    df = pd.read_csv(per_read_modifications, delimiter='\t')

    # Map each field to its own column using the names in the dictionary
    df['reference_sequence'] = df['read_id'].map(lambda read_id: primary_alignments.get(read_id, {}).get('reference_sequence'))
    df['flag'] = df['read_id'].map(lambda read_id: primary_alignments.get(read_id, {}).get('flag'))
    df['ref_start_pos'] = df['read_id'].map(lambda read_id: primary_alignments.get(read_id, {}).get('ref_start_pos'))
    df['cigar_string'] = df['read_id'].map(lambda read_id: primary_alignments.get(read_id, {}).get('cigar_string'))
    df['ref_length'] = df['read_id'].map(lambda read_id: primary_alignments.get(read_id, {}).get('reference_length'))
    df['basecalled_sequence'] = df['read_id'].map(lambda read_id: primary_alignments.get(read_id, {}).get('basecalled_sequence'))
    df['q_score'] = df['read_id'].map(lambda read_id: primary_alignments.get(read_id, {}).get('q_score'))

   # Write the updated DataFrame back to a new TSV file
    df.to_csv(per_read_modifications, sep='\t', index=False)

    return per_read_modifications

def process_per_read_mod_results(mod_base, can_base, analysis_dir, per_read_mod, min_coverage=100):
    """
    Processes the per-read modification results by grouping data based on the reference sequence, 
    creating subdirectories for each reference, and calling the data_analysis function on each subset.

    Parameters:
        mod_base (str): Modified base used in the analysis.
        can_base (str): Canonical base used in the analysis.
        analysis_dir (str): Path to the directory where the analysis will be stored.
        per_read_mod (str): Path to the per-read modification TSV file.
        min_coverage (int): Minimum number of unique reads required for processing (default is 100).

    Returns:
        None
    """
    # Read the TSV file into a DataFrame
    df = pd.read_csv(per_read_mod, sep='\t')

    # Create the metadata DataFrame with two columns: 'XNADNA_model' and 'meta_info'
    df_ref_pass = pd.DataFrame(columns=["passed_reference_sequences"])
    
    # Get unique reference sequences
    unique_references = df['reference_sequence'].unique()

    print(f'Xemora [STATUS] - Performing consensus level analysis on reference sequences with more than {min_coverage} read coverage.')
    for reference in unique_references:
        # Filter rows where reference_sequence column == reference
        subset_df = df[df['reference_sequence'] == reference]

        # Calculate alignment coverage
        alignment_coverage = subset_df['read_id'].nunique() 

        if alignment_coverage >= min_coverage:
            # Report on alignment and coverage
            print(f'Xemora [STATUS] - Analyzing sequences mapping to reference = {reference} | coverage {alignment_coverage}')
            
            # Append reference-specific metadata in the 'passed_reference_sequences column
            df_ref_pass = pd.concat([df_ref_pass, pd.DataFrame({"passed_reference_sequences": [f'ref_{reference}']})], ignore_index=True)

            # Create a directory using the actual reference sequence name
            dir_name = os.path.join(analysis_dir, f'ref_{reference}')
            os.makedirs(dir_name, exist_ok=True)

            # Save the subset DataFrame into the subdirectory
            output_file_path = os.path.join(dir_name, f'per_read_summary_{reference}.tsv')
            subset_df.to_csv(output_file_path, sep='\t', index=False)

            # Call the data analysis function
            data_analysis(dir_name, output_file_path)
        else:
            continue

        # Save df_meta as metadata in the analysis directory
        passed_reference_files = os.path.join(analysis_dir, 'passed_references.tsv')
        df_ref_pass.to_csv(passed_reference_files, sep='\t', index=False)

def data_analysis(analysis_dir, per_read_modifications):
    """
    This function performs two key analyses on the provided per-read modification data.

    Step 1: Realignment/Normalization
    - Uses the script 'xr_realign.py' to realign the per-read modification data.
    - The output is saved as 'per_read_modifications_normalized.tsv' in the specified analysis directory.

    Step 2: Focus Position Analysis and Plot Generation
    - Uses the script 'xr_focus_position.py' to perform focus base modification analysis.
    - Two output files are generated:
      1. A plot showing combined base modification data.
      2. A plot focusing on odds of modification occurrences.

    Parameters:
        analysis_dir (str): The directory where the analysis results will be saved.
        per_read_modifications (str): Path to the per-read modification TSV file.

    Returns:
        None
    """
    normalized_per_read_modifications = os.path.join(analysis_dir, 'per_read_modifications_normalized.tsv')
    cmd = 'python lib/xr_realign.py '+per_read_modifications+' '+normalized_per_read_modifications
    os.system(cmd)

    focus_normalized_per_read_modifications = os.path.join(os.path.dirname(normalized_per_read_modifications), 'focus_base_modifications_normalized.tsv')
    graph_path = os.path.join(analysis_dir, 'combined_plots.pdf')
    odds_graph_path = os.path.join(analysis_dir, 'odds_plots.pdf')
    entropy_graph_path = os.path.join(analysis_dir, 'entropy_plots.pdf')
    cmd = 'python lib/xr_focus_position.py '+normalized_per_read_modifications+' '+focus_normalized_per_read_modifications+' '+graph_path+' '+odds_graph_path+' '+entropy_graph_path
    os.system(cmd)
    return focus_normalized_per_read_modifications

def plot_overlay(focus_base_list):
    """
    takes in a list of file pathways containing the analyzed focus position data

    WIP
    """
    return focus_base_list

def check_consensus_directories(analysis_dir):
    """
    Verifies that the appropriate model directories and metadata files exist in the analysis directory.

    Parameters:
        analysis_dir (str): The path to the base analysis directory.
        
    Returns:
        bool: True if all directories and metadata files are present, False otherwise.
    """

    # Define the expected directory names
    model_dirs = {
        'model_1_fwd': os.path.join(analysis_dir, 'model_1_fwd'),
        'model_1_rev': os.path.join(analysis_dir, 'model_1_rev'),
        'model_2_fwd': os.path.join(analysis_dir, 'model_2_fwd'),
        'model_2_rev': os.path.join(analysis_dir, 'model_2_rev')
    }

    # Check for the presence of each directory and its metadata.tsv file
    for model_name, model_dir in model_dirs.items():
        if not os.path.isdir(model_dir):
            print(f'Xemora [ERROR] - Directory missing: {model_dir}')
            return False

        # Check for the presence of the metadata.tsv file
        passed_ref_file = os.path.join(model_dir, 'passed_references.tsv')
        if not os.path.exists(passed_ref_file):
            print(f'Xemora [ERROR] - Passed references file missing in {model_dir}')
            return False

    print('Xemora [STATUS] - All required model directories and file containing passed reference sequences are present.')
    return True

def load_passed_references(passed_ref_dir):
    """
    Loads the passed references from the passed_references.tsv file.

    Parameters:
        passed_ref_dir (str): Directory where the passed_references.tsv file is stored.
        
    Returns:
        set: A set of reference sequences that passed the filtering criteria.
    """
    passed_ref_file = os.path.join(passed_ref_dir, 'passed_references.tsv')
    if os.path.exists(passed_ref_file):
        df = pd.read_csv(passed_ref_file, sep='\t')
        return set(df['passed_reference_sequences'])  # Returning a set for fast lookup
    else:
        raise FileNotFoundError(f'Xemora [ERROR] - Passed reference file missing in {passed_ref_dir}')


def pair_forward_reverse(analysis_dir):
    """
    Pairs forward strand data from model_1_fwd with reverse strand data from model_2_rev,
    and reverse strand data from model_1_rev with forward strand data from model_2_fwd.
    
    The pairing is done based on matching reference sequences that passed the filtering criteria.
    
    Parameters:
        analysis_dir (str): Path to the base analysis directory containing model directories.

    Returns:
        dict: A dictionary where the keys are reference sequences and the values are dictionaries with
              file paths for forward and reverse strands, in the form:
              {
                  "fwd_rev": {
                      "fwd": path_to_model_1_fwd_file,
                      "rev": path_to_model_2_rev_file
                  },
                  "rev_fwd": {
                      "fwd": path_to_model_2_fwd_file,
                      "rev": path_to_model_1_rev_file
                  }
              }
    """
    # Define model directories
    model_dirs = {
        'model_1_fwd': os.path.join(analysis_dir, 'model_1_fwd'),
        'model_1_rev': os.path.join(analysis_dir, 'model_1_rev'),
        'model_2_fwd': os.path.join(analysis_dir, 'model_2_fwd'),
        'model_2_rev': os.path.join(analysis_dir, 'model_2_rev')
    }

    # Load passed reference sequences from both model_1_fwd and model_2_rev
    passed_references_model_1_fwd = load_passed_references(model_dirs['model_1_fwd'])
    passed_references_model_2_rev = load_passed_references(model_dirs['model_2_rev'])
    passed_references_model_1_rev = load_passed_references(model_dirs['model_1_rev'])
    passed_references_model_2_fwd = load_passed_references(model_dirs['model_2_fwd'])

    # Combine all passed references
    passed_references = passed_references_model_1_fwd.intersection(passed_references_model_2_rev).union(
        passed_references_model_1_rev.intersection(passed_references_model_2_fwd)
    )

    # Initialize a dictionary to store the forward/reverse pairings
    paired_results = {}

    # Iterate over passed references and find pairings for forward and reverse strands
    for ref in passed_references:
        ref_group = {}

        # Construct paths for model_1_fwd and model_2_rev (fwd vs rev)
        fwd_model_1 = os.path.join(model_dirs['model_1_fwd'], ref, 'focus_base_modifications_normalized.tsv')
        rev_model_2 = os.path.join(model_dirs['model_2_rev'], ref, 'focus_base_modifications_normalized.tsv')

        if os.path.exists(fwd_model_1) and os.path.exists(rev_model_2):
            ref_group["fwd_rev"] = {
                "fwd": fwd_model_1,
                "rev": rev_model_2
            }

        # Construct paths for model_1_rev and model_2_fwd (rev vs fwd)
        rev_model_1 = os.path.join(model_dirs['model_1_rev'], ref, 'focus_base_modifications_normalized.tsv')
        fwd_model_2 = os.path.join(model_dirs['model_2_fwd'], ref, 'focus_base_modifications_normalized.tsv')

        if os.path.exists(rev_model_1) and os.path.exists(fwd_model_2):
            ref_group["rev_fwd"] = {
                "fwd": fwd_model_2,
                "rev": rev_model_1
            }

        # Only add to the results if both pairings are found
        if ref_group:
            paired_results[ref] = ref_group

    return paired_results
    
def pairwise_analysis(paired_results, pairwise_dir):
    """
    Pairwise analysis function to iterate through the paired_results dictionary,
    generate directories for each reference, and run the pairwise statistics command.

    Parameters:
        paired_results (dict): The dictionary containing paired forward and reverse strand file paths.
        base_output_dir (str): The base directory where all pairwise analysis results will be saved.
    
    Returns:
        None
    """

    # Step 2: Iterate over the paired results dictionary
    for ref, pairings in paired_results.items():
        # Step 3: Create a per-reference directory
        ref_output_dir = os.path.join(pairwise_dir, ref)
        
        # Step 4: Handle the fwd_rev pairing (model_1_fwd:model_2_rev)
        if "fwd_rev" in pairings:
            fwd_rev_dir = os.path.join(ref_output_dir, 'model_1_fwd_model_2_rev')

            model_1_statistics = pairings["fwd_rev"]["fwd"]
            model_2_statistics = pairings["fwd_rev"]["rev"]

            # Construct and run the command
            cmd = f'python lib/xr_pairwise_stats.py {fwd_rev_dir} {model_1_statistics} {model_2_statistics}'
            print(cmd)
            os.system(cmd)

        # Step 5: Handle the rev_fwd pairing (model_2_fwd:model_1_rev)
        if "rev_fwd" in pairings:
            rev_fwd_dir = os.path.join(ref_output_dir, 'model_2_fwd_model_1_rev')

            model_1_statistics = pairings["rev_fwd"]["fwd"]
            model_2_statistics = pairings["rev_fwd"]["rev"]

            # Construct and run the command
            cmd = f'python lib/xr_pairwise_stats.py {rev_fwd_dir} {model_1_statistics} {model_2_statistics}'
            print(cmd)
            os.system(cmd)

    print('Xemora [STATUS] - Pairwise analysis completed.')
    
def main():
    """
    main runs methods to perform XNA basecalling. 
    """

    #Set paramters if present
    set_params(params)

    #Check FASTA file for XNA. If found, perform hard swap (method 1). If none found, perform full read inference (method 2) 
    xna_present = check_fasta(ref_fasta, xna_base_pairs) # xna_base_pairs present in xr_params

    if xna_present:
        # Method 1 - XNA hardswap
        # xFASTA generation
        fasta_file = xfasta_conversion(ref_fasta, ref_dir)
    elif not xna_present:
        # Method 2 - Full reference analysis
        fasta_file, ref_lengths = copy_fasta(ref_fasta, ref_dir)
    else:
        # Future work placeholder for de novo basecalling
        print('Xemora [STATUS] - Placeholder for method 3, de novo basecalling. Please check file inputs.')
        sys.exit()

    #Verify only one file type in raw data directory
    filetype = validate_read_directory(raw_dir)
    
    #Perform either FAST5 converison or POD55 merging
    if filetype == 'fast5':
        merged_pod5 = cod5_to_fast5(raw_dir, pod_dir, overwrite_pod)
    elif filetype == 'pod5':
        merged_pod5 = pod5_merge(raw_dir, pod_dir, overwrite_pod)
    else:
        print('Xemora [ERROR] - file type in raw data directory is not POD5 or FAST5, please check if your directory is correct') 
        sys.exit()
    
    #Perform basecalling using Dorado
    bc_bam = dorado_basecall(os.path.expanduser(dorado_path), dorado_model, min_qscore, merged_pod5, mod_bam_dir, basecall_pod, max_reads, filter_readIDs)
    
    #Extra filtering on basecalls
    print(f'Xemora [STATUS] - min_qscore_analysis : {min_qscore_analysis } - Performing filtering of reads by q-score')
    bc_bam = filter_bam_by_qs(bc_bam, mod_bam_dir)

    #Align using minimap2 
    aligned_bam = minimap2_aligner(bc_bam, fasta_file, mod_bam_dir)

    #AS filter 
    filter_by_alignment_score(aligned_bam)

    #Get primary alignments for reference sequence assignment
    primary_alignments = get_primary_alignments(aligned_bam)

    #Initialize focus position file path
    focus_position_list = []


    #Loop through model files and perform Remora validate/ data analysis
    for i in range(len(model_files)):
        directions = ['fwd', 'rev']
        model_file = model_files[i]

        #Non-zero naming scheme for user ease of understanding 
        model_number = i + 1 
    #for model_file in model_files: 
        mod_base, can_base = load_xna_model(model_file)
        
        # Create separate directories for forward and reverse strands
        '''
            maybe something like put xna and model number in bed name X_model_1_fwd.bed or something 
        ''' 
        #Bed file generation logic
        if xna_present:
            # Method 1 - XNA hardswap
            # Generate BED files for regions to analyze (both forward and reverse)
            bed_list = bed_gen(fasta_file, mod_base, mod_base, mod_chunk_range, mod_chunk_shift, model_number)
            
        elif not xna_present:
            # Method 2 - Full reference analysis

            # Generate BED files for the full reference for both forward and reverse
            bed_list = full_read_bed(ref_lengths, mod_base, mod_base, ref_dir)
            
        else:
            # Future work placeholder for de novo basecalling
            print('Xemora [STATUS] - Placeholder for method 3, de novo basecalling. Please check file inputs.')
            sys.exit()

        # Process forward and reverse files together
        for direction, bed_file in zip(directions, bed_list):
            # Generate chunk files to analyze 
            basecalling_chunks = generate_chunks(merged_pod5, aligned_bam, os.path.join(chunk_dir,f'basecalling_chunks_{model_number}_{direction}.npz'), bed_file, mod_base, kmer_context, kmer_table_path, regenerate_chunks) #changed to take in a chunk file

            if run_analysis or not os.path.exists(os.path.join(remora_output_dir, f'per-read_modifications_model_{model_number}_{direction}.tsv')) or not os.path.exists(os.path.join(remora_output_dir,f'summary_modifications_model_{model_number}_{direction}.tsv')):
                # Run Xemora basecall for this chunk
                per_read_mod, summary_mod = xemora_basecall(basecalling_chunks, model_file, os.path.join(remora_output_dir,f'per-read_modifications_model_{model_number}_{direction}.tsv'), os.path.join(remora_output_dir,f'summary_modifications_model_{model_number}_{direction}.tsv'))
                print(per_read_mod)

                # Add reference sequence to per-read modifications
                per_read_mod = add_per_read_mapping(primary_alignments, per_read_mod)
                print(per_read_mod)

                # Perform analysis for each individual sequence
                '''
                problems start happening here, need a better naming scheme than 'ref_1, ref_2, etc'
                here is where we might want to add a sub model directory (e.g., analysis/model_1/ref_1)
                '''
                per_model_analysis_dir = check_make_dir(os.path.join(analysis_dir, f'model_{model_number}_{direction}'))
                focus_position_results = process_per_read_mod_results(mod_base, can_base, per_model_analysis_dir, per_read_mod)  # `min_cluster_size` present in xr_params
                focus_position_list.append(focus_position_results)
            else: 
                print('Xemora [STATUS] - Skipping Remora validate and consensus sequence analysis')

    if len(model_files) == 2 and pairwise: 
        print('Xemora [STATUS] - Performing pairwise analysis')
        '''
        need: function that does the association (B fwd, S rev; S fwd, B rev)
        metadata 
        actually running the pairwise stats script 

        need to add stopping logic if more than 4 directories present. need metadata files to know what models were.
        '''
        check_con_dir = check_consensus_directories(analysis_dir)
        
        if check_con_dir:
            #Pair reference sequences for corresponding models 
            ref_pairings = pair_forward_reverse(analysis_dir)
            print(ref_pairings)

            #Run pairwise analysis: 
            pairwise_dir = check_make_dir(os.path.join(working_dir, 'pairwise_results'))
            pairwise_analysis(ref_pairings, pairwise_dir)
        else:
            print('Xemora [ERROR] - Passed references file not found')
    elif not pairwise:
        print('Xemora [STATUS] - Pairwise flag set to False, skipping pairwise analysis')
    elif len(model_files) > 2 or len(model_files) < 2 and pairwise:
        print(f'Xemora [ERROR] - Tried to perform pairwise analysis with {len(model_files)}. Can only perform pairwise analysis with 2 model file inputs')
        sys.exit()
    else:
        print('Xemora [STATUS] - Skipping pairwise analysis')

    print('Xemora Basecall [COMPLETE] - Have a nice day :)')

if __name__ == "__main__":
    main()
