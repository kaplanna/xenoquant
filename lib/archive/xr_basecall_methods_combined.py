########################################################################
########################################################################
"""
xr_basecall_methods.py 

Title: Unpublished work

Perform 'standard' single context basecalling as done prior except this version of the 
script uses dorado and has the start for the single context consensus analysis 
This script is method based and the main version of the single context / 
reference anchored basecall

Note: if you need to troubleshoot/change the consensus analysis, the script is 

'xr_consensus_modifications.py' and the output file to analyze in the working directory output_modification_results.tsv

By: H. Kawabe, N. Kaplan, J. Sumabat, J. A. Marchand

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
import pysam
import pandas as pd

############################################################
working_dir = os.path.expanduser(sys.argv[1])
raw_dir = os.path.expanduser(sys.argv[2])
xna_ref_fasta = os.path.expanduser(sys.argv[3])
model_file = os.path.expanduser(sys.argv[4])

#########################################################################
#Generate directories
print('Xemora [Status] - Initializing Xemora basecalling.')

working_dir = check_make_dir(working_dir)
ref_dir = check_make_dir(os.path.join(working_dir,'references'))
chunk_dir = check_make_dir(os.path.join(working_dir,'chunks'))
mod_dir = check_make_dir(os.path.join(working_dir,'preprocess'))
mod_pod_dir = check_make_dir(os.path.join(mod_dir,'pod5'))
mod_bam_dir = check_make_dir(os.path.join(mod_dir,'bam'))
remora_output_dir = check_make_dir(os.path.join(working_dir, 'remora_outputs'))
#########################################################################

def get_revcomp_base(base, base_pairs):
    """
    Given a base and a list of 2-letter base pairs, return the complementary base.
    """
    for pair in base_pairs:
        if base in pair:
            return pair.replace(base, "")
    raise ValueError(f"Complement for base '{base}' not found in provided base pairs.")

#Methods 

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
        # Add Dorado summary step here after basecalling
        summary_cmd = 'dorado summary -v {} > {}'.format(output_bam, os.path.join(bam_directory, 'sequencing_summary.txt'))
        os.system(summary_cmd)
        print('Xemora [STATUS] - Dorado summary saved to sequencing_summary.txt.')
        
        return output_bam
        
    else:
        print('Xemora [STATUS] - Skipping POD5 basecalling for modified bases.')
        return output_bam


def minimap2_aligner(input_bam, xfasta_path, bam_directory):
    """
    minimap2_aligner takes in a bam file from Dorado basecaller and outputs a new 
    aligned bam file using an xFASTA formatted file 

    Parameters:
    input_bam - file path to bam file to align
    xfasta_path - file pathway to xFASTA file to align basecalls to 
    bam_directory - output file pathway for aligned BAM file 

    Returns: 
    aligned_bam - minimap2 aligned bam file 
    """
    aligned_bam = os.path.join(bam_directory, 'aligned.BAM')
    cmd = 'samtools fastq -T "*" '+input_bam+ ' | minimap2 -y -ax map-ont --score-N 0 --secondary no --sam-hit-only --MD '+xfasta_path+ ' - | samtools view -F0x800 -bho ' +aligned_bam
    os.system(cmd)
    
    return aligned_bam

    return aligned_bam
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
#Step 4: Bed file generation 
def bed_gen(input_fasta, xna_base, sub_base, chunk_range, chunk_shift): 
    """
    need to implement shift (maybe)
    
    bed_gen takes in an xFASTA file, xna base, a base to substitute with the 
    xna base, region for analysis around focus, and shift from focus to generate 
    a BED file with the region to perform Remora chunk generation on. 
    
    Parameters:
    input_fasta - xFASTA file path 
    xna_base - xna of interests, in xr_params 
    sub_base - base to substitute
    chunk_range - integer +/- around focus to analyze 
    chunk_shift - integer with how many bases to move the focus region
    
    Returns:
    output_bed - BED file pathway that is generated 
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

def generate_chunks(pod_file, bam_file, chunk_dir, bed_file, mod_base, kmer_context, kmer_table_path, regenerate_chunks):
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
    """
    #Step 5: Generate Chunks. 
    chunk_file = os.path.join(chunk_dir, 'basecall_chunks.npz')
    if regenerate_chunks == True:
        '''
        print('Xemora  [STATUS] - Generating chunks for modified basecalling.')
        cmd = 'remora \
          dataset prepare \
          '+os.path.join(mod_pod_dir,os.path.basename(raw_dir))+'.pod5'+' \
          '+mod_bam_path+'_all.bam'+' \
          --output-remora-training-file '+os.path.join(chunk_dir,'basecall_chunks.npz')+' \
          --focus-reference-positions '+os.path.join(ref_dir,mod_base)+'.bed'+' \
          --mod-base '+mod_base+' '+mod_base+' \
          --motif '+can_base+' 0 \
          --kmer-context-bases '+kmer_context+' \
          --refine-kmer-level-table '+kmer_table_path+' \
          --refine-rough-rescale '+' \
          --chunk-context '+chunk_context
        os.system(cmd)
        '''
        print('Xemora  [STATUS] - Generating chunks for modified basecalling.')
        #no motif
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
        print('Xemora [STATUS] - Skipping chunk generation')
        return chunk_file

def xemora_basecall(working_dir, chunk_file, model_file):
    """
    xemora_basecall takes in the generate chunk file and a Xemora model and performs 
    inference on the chunk file to determine is XNAs are present 
    
    Parameters: 
    working_dir - output directory inference files
    chunk_file - file pathway to a Remora chunk file 
    model_file - file pathway to a Xemora model 
    
    Returns: 
    """
    print('Xemora  [STATUS] - Performing basecalling.')
    try:
        per_read_output = os.path.join(remora_output_dir, 'per-read_modifications.tsv')
        summary_output = os.path.join(remora_output_dir, 'summary_modifications.tsv')
        
        cmd = 'remora \
          validate from_remora_dataset \
          '+chunk_file+' \
          --model '+os.path.expanduser(model_file)+' \
          --full-results-filename '+per_read_output+' \
          --out-file '+summary_output
        os.system(cmd)

        print('Xemora [STATUS] - Basecalling done.')
        print('Xemora [STATUS] - Basecalling done. Saving results '+per_read_output)
        print('Xemora [STATUS] - Basecalling done. Saving results '+summary_output)
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

# Updated main()
def main():
    from xr_params import mod_base as top_mod_base, can_base as top_can_base


    # Shared steps (run once)
    xfasta_file = xfasta_conversion(xna_ref_fasta, ref_dir)
    filetype = validate_read_directory(raw_dir)

    if filetype == 'fast5':
        merged_pod5 = cod5_to_fast5(raw_dir, mod_pod_dir, overwrite_pod)
    elif filetype == 'pod5':
        merged_pod5 = pod5_merge(raw_dir, mod_pod_dir, overwrite_pod)
    else:
        print('Xemora [ERROR] - Invalid raw data format.')
        sys.exit()

    bc_bam = dorado_basecall(dorado_path, dorado_model, min_qscore, merged_pod5, mod_bam_dir, basecall_pod, max_bc_reads, filter_readIDs_bc)
    aligned_bam = minimap2_aligner(bc_bam, xfasta_file, mod_bam_dir)
    primary_alignments = get_primary_alignments(aligned_bam)

    # Strand logic: derive bottom strand mod_base
    bottom_mod_base = get_revcomp_base(top_mod_base, xna_base_pairs)
    bottom_can_base = get_revcomp_base(top_can_base, standard_base_pairs)

    strand_configs = [
        {"mod_base": top_mod_base, "can_base": top_can_base},
        {"mod_base": bottom_mod_base, "can_base": bottom_can_base}
    ]

    for config in strand_configs:
        mod_base = config["mod_base"]
        can_base = config["can_base"]
        print(f'\nXemora [STATUS] - Running mod_base: {mod_base}')

        strand_bed_file = bed_gen(
            xfasta_file,
            mod_base,
            can_base,
            mod_chunk_range,
            mod_chunk_shift
        )

        strand_chunk_file = generate_chunks(
            merged_pod5,
            aligned_bam,
            chunk_dir,
            strand_bed_file,
            mod_base,
            kmer_context,
            kmer_table_path,
            regenerate_chunks
        )

        per_read_mod, summary_mod = xemora_basecall(working_dir, strand_chunk_file, model_file)
        per_read_mod = add_per_read_mapping(primary_alignments, per_read_mod)

        # Save results using mod_base in the filename
        mod_base_suffix = f'_{mod_base}'
        per_read_mod_out = per_read_mod.replace('.tsv', f'{mod_base_suffix}.tsv')
        summary_mod_out = summary_mod.replace('.tsv', f'{mod_base_suffix}.tsv')
        os.rename(per_read_mod, per_read_mod_out)
        os.rename(summary_mod, summary_mod_out)

        print(f'Xemora [STATUS] - Saved results for {mod_base}')

if __name__ == "__main__":
    main()

