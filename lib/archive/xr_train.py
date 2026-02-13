########################################################################
########################################################################
"""
xr_train.py 

Title: Unpublished work

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


############################################################
print('Xemora [Status] - Initializing Xemora...')

#Initialize
working_dir = os.path.expanduser(sys.argv[1])
xna_raw_dir = os.path.expanduser(sys.argv[2])
xna_ref_fasta = os.path.expanduser(sys.argv[3])
dna_raw_dir = os.path.expanduser(sys.argv[4])
dna_ref_fasta = os.path.expanduser(sys.argv[5])
#ext_dir = "~/xenomorph/xemora-beta/xemora-test/230309_Xemora_PZ_GC_train/chunks/training_chunks.npz"

#Generate directories
working_dir = check_make_dir(working_dir)
ref_dir = check_make_dir(os.path.join(working_dir,'references'))
model_dir = check_make_dir(os.path.join(working_dir,'model'))
chunk_dir = check_make_dir(os.path.join(working_dir,'chunks'))
mod_dir = check_make_dir(os.path.join(working_dir,'modified'))
mod_pod_dir = check_make_dir(os.path.join(mod_dir,'pod5'))
mod_fastq_dir = check_make_dir(os.path.join(mod_dir,'fastq'))
mod_bam_dir = check_make_dir(os.path.join(mod_dir,'bam'))
can_dir = check_make_dir(os.path.join(working_dir,'canonical'))
can_pod_dir = check_make_dir(os.path.join(can_dir,'pod5'))
can_fastq_dir = check_make_dir(os.path.join(can_dir,'fastq'))
can_bam_dir = check_make_dir(os.path.join(can_dir,'bam'))
############################################################


#Step 0: FASTA to xFASTA conversion
if os.path.isfile(os.path.expanduser(xna_ref_fasta)): 
    cmd = 'python lib/xr_fasta2x_rc.py '+os.path.expanduser(xna_ref_fasta)+' '+os.path.join(ref_dir,'x'+os.path.basename(xna_ref_fasta))
    os.system(cmd)
else: 
    print('Xemora [ERROR] - XNA Reference fasta file not found. Please check file exist or file path.')
    sys.exit()

if os.path.isfile(os.path.expanduser(dna_ref_fasta)): 
    cmd = 'python lib/xr_fasta2x_rc.py '+os.path.expanduser(dna_ref_fasta)+' '+os.path.join(ref_dir,'x'+os.path.basename(dna_ref_fasta))
    os.system(cmd)
else: 
    print('Xemora [ERROR] - DNA Reference fasta file not found. Please check file exist or file path.')
    sys.exit()

#Step 1: Generate or merge pod5 files if needed
xna_file_type = os.listdir(xna_raw_dir)
# Check if the directory is not empty
if xna_file_type:
    # Get the first file in the directory
    xna_first_file = xna_file_type[0]
    # Check if the first file is a .pod5 file
    if xna_first_file.endswith(".pod5"):
        if os.path.isfile(os.path.join(mod_pod_dir,os.path.basename(xna_raw_dir))+'.pod5')==False:
            pod5_merge(get_pod5_subdir(xna_raw_dir), os.path.join(mod_pod_dir,os.path.basename(xna_raw_dir))+'.pod5')
        else:
            print('Xemora [STATUS] - POD5 files merged. Skipping merge')
    else:
        if os.path.isfile(os.path.join(mod_pod_dir,os.path.basename(xna_raw_dir))+'.pod5')==False:
            cod5_to_fast5(get_fast5_subdir(xna_raw_dir), os.path.join(mod_pod_dir,os.path.basename(xna_raw_dir))+'.pod5')
        else: 
            print('Xemora [STATUS] - Converted POD5 file for modified base found. Skipping POD5 coversion')
else:
    print('Xemora [ERROR] - Modified Fast5/POD5 directory empty, please check input directory')
    sys.exit()

dna_file_type = os.listdir(dna_raw_dir)
# Check if the directory is not empty 
if dna_file_type:
    # Get the first file in the directory 
    dna_first_file = dna_file_type[0]
    #check if the first file is a pod5 file
    if dna_first_file.endswith(".pod5"):
        if os.path.isfile(os.path.join(can_pod_dir,os.path.basename(dna_raw_dir))+'.pod5')==False:
            pod5_merge(get_pod5_subdir(dna_raw_dir), os.path.join(can_pod_dir,os.path.basename(dna_raw_dir))+'.pod5')
        else:
            print('Xemora [STATUS] - POD5 files merged. Skipping merge')
    else: 
        if os.path.isfile(os.path.join(can_pod_dir,os.path.basename(dna_raw_dir))+'.pod5')==False: 
            cod5_to_fast5(get_fast5_subdir(dna_raw_dir), os.path.join(can_pod_dir,os.path.basename(dna_raw_dir))+'.pod5')
        else: 
            print('Xemora [STATUS] - Converted POD5 file for canonical base found. Skipping POD5 coversion')
else:
    print('Xemora [ERROR] - Canonical Fast5/POD5 directory empty, please check input directory')
    sys.exit()
    
#Step 2: #Basecall pod5 files 
if basecall_pod ==True: 
    cmd=os.path.expanduser(basecaller_path)+' -i '+mod_pod_dir+' -s '+mod_fastq_dir+' -c '+guppy_config_file+' -x auto --bam_out --index --moves_out -a '+os.path.join(ref_dir,'x'+os.path.basename(xna_ref_fasta))+' --min_qscore '+str(training_mod_min_qscore)
    os.system(cmd)
else: 
    print('Xemora [STATUS] - Skipping POD5 basecalling for modified bases.')


if basecall_pod ==True:
    cmd=os.path.expanduser(basecaller_path)+' -i '+can_pod_dir+' -s '+can_fastq_dir+' -c '+guppy_config_file+' -x auto --bam_out --index --moves_out -a '+os.path.join(ref_dir,'x'+os.path.basename(dna_ref_fasta))+' --min_qscore '+str(training_can_min_qscore)
    os.system(cmd)
else: 
    print('Xemora [STATUS] - Skipping POD5 basecalling for canonical bases.')


################################################################################################################################

# Merge modified bam files
if os.path.isfile(os.path.join(mod_bam_dir, os.path.basename(mod_bam_dir)) + '.bam') == False or regenerate_bam == True: 

    #Purge directory
    cmd_rmv = 'rm -rf '+mod_bam_dir+'/*'
    os.system(cmd_rmv)
    
    #Set modified bam directory path 
    mod_bam_path = os.path.join(mod_bam_dir, os.path.basename(mod_bam_dir))
    
    # Merging pass bam files
    cmd_pass = 'samtools merge ' + mod_bam_path + '_pass.bam ' + os.path.join(mod_fastq_dir, 'pass/*.bam -f')
    print('Xemora [STATUS] - Merging modified PASS BAM files.')
    os.system(cmd_pass)
    
    print('Amount of total reads in modified PASS BAM file') 
    cmd = 'samtools view -c ' + mod_bam_path + '_pass.bam'
    os.system(cmd) 
    print('Amount of unaligned reads in modified PASS BAM file') 
    cmd = 'samtools view -c -f 4 ' + mod_bam_path + '_pass.bam'
    os.system(cmd)
    
    
    if merge_fail == True:
        # Merging fail bam files
        print('Xemora [STATUS] - Merging modified FAIL BAM files.')
        cmd_fail = 'samtools merge ' + mod_bam_path + '_fail.bam ' + os.path.join(mod_fastq_dir, 'fail/*.bam -f')
        os.system(cmd_fail)
        
        print('Amount of total reads in modified FAIL BAM file')
        cmd = 'samtools view -c ' + mod_bam_path + '_fail.bam'
        os.system(cmd) 
        print('Amount of unaligned reads in modified FAIL BAM file')
        cmd = 'samtools view -c -f 4 ' + mod_bam_path + '_fail.bam'
        os.system(cmd)
        
        
    #Generate merged bam
    cmd_both = 'samtools merge ' + mod_bam_path + '_all.bam ' + mod_bam_path+'*.bam -f'
    os.system(cmd_both)

    print('Amount of total reads in modified FULL BAM file') 
    cmd = 'samtools view -c ' + mod_bam_path + '_all.bam'
    os.system(cmd) 
    print('Amount of unaligned reads in modified FULL BAM file') 
    cmd = 'samtools view -c -f 4 ' + mod_bam_path + '_all.bam'
    os.system(cmd)



# Merge canonical bam files
if os.path.isfile(os.path.join(can_bam_dir, os.path.basename(can_bam_dir)) + '.bam') == False or regenerate_bam == True: 

    #Clear directory
    cmd_rmv = 'rm -rf '+can_bam_dir+'/*'
    os.system(cmd_rmv)

    #Set canonical bam directory path 
    can_bam_path = os.path.join(can_bam_dir, os.path.basename(can_bam_dir))
    
    #Merging pass bam files
    cmd_pass = 'samtools merge ' + can_bam_path + '_pass.bam ' + os.path.join(can_fastq_dir, 'pass/*.bam -f')
    print('Xemora [STATUS] - Merging canonical PASS BAM files.')
    os.system(cmd_pass)

    print('Xemora [STATUS] - Indexing canonical PASS BAM file') 
    cmd_index = 'samtools index ' + can_bam_path + '_pass.bam'
    os.system(cmd_index)
    
    print('Amount of total reads in canonical PASS BAM file') 
    cmd = 'samtools view -c ' + can_bam_path + '_pass.bam'
    os.system(cmd) 
    print('Amount of unaligned reads in canonical PASS BAM file') 
    cmd = 'samtools view -c -f 4 ' + can_bam_path + '_pass.bam'
    os.system(cmd)

    if merge_fail == True:

        # Merging fail bam files
        print('Xemora [STATUS] - Merging canonical FAIL BAM files.')
        cmd_fail = 'samtools merge ' + can_bam_path + '_fail.bam ' + os.path.join(can_fastq_dir, 'fail/*.bam -f')
        os.system(cmd_fail)

        print('Xemora [STATUS] - Indexing canonical FAIL BAM file') 
        cmd_index = 'samtools index ' + can_bam_path + '_fail.bam'
        os.system(cmd_index)
        
        print('Amount of total reads in canonical FAIL BAM file') 
        cmd = 'samtools view -c ' + can_bam_path + '_fail.bam'
        os.system(cmd) 
        print('Amount of unaligned reads in canonical FAIL BAM file') 
        cmd = 'samtools view -c -f 4 ' + can_bam_path + '_fail.bam'
        os.system(cmd)
        
    #Generate merged bam
    cmd_both = 'samtools merge ' + can_bam_path + '_all.bam ' + can_bam_path+'*.bam -f'
    os.system(cmd_both)
    
    print('Xemora [STATUS] - Indexing canonical full BAM file') 
    cmd_index = 'samtools index ' + can_bam_path + '_all.bam'
    os.system(cmd_index)
    
    print('Amount of total reads in canonical FULL BAM file') 
    cmd = 'samtools view -c ' + can_bam_path + '_all.bam'
    os.system(cmd) 
    print('Amount of unaligned reads in canonical FULL BAM file') 
    cmd = 'samtools view -c -f 4 ' + can_bam_path + '_all.bam'
    os.system(cmd)
    
    #BAM primary alignment only filtering 
    print('Xemora [STATUS] - Removing unaligned, secondary, and supplementary alignments')
    filter_primary_alignments(mod_bam_path + '_all.bam', mod_bam_path + '_filtered.bam') #modified dataset filtering 
    filter_primary_alignments(can_bam_path + '_all.bam', can_bam_path + '_filtered.bam')
################################################################################################################################




################################################################################################################################
#Step 4: Bed file generation 
if os.stat(os.path.join(ref_dir,'x'+os.path.basename(xna_ref_fasta))).st_size == 0: 
    print('Xemora [ERROR] - Empty xfasta file generated. Check that XNA bases were present in sequence of input fasta file.')
    sys.exit()

print('Xemora [STATUS] - Generating bed file for modified base.')
cmd = 'python lib/xr_xfasta2bed.py '+os.path.join(ref_dir,'x'+os.path.basename(xna_ref_fasta))+' '+os.path.join(ref_dir,mod_base+'.bed ' +mod_base+' '+mod_base)
os.system(cmd)

print('Xemora [STATUS] - Generating bed file for canonical base.')
cmd = 'python lib/xr_xfasta2bed.py '+os.path.join(ref_dir,'x'+os.path.basename(dna_ref_fasta))+' '+os.path.join(ref_dir,can_base+'.bed '+mod_base+' '+can_base)
os.system(cmd)

#Optional Bed Filtering
if bed_filtering == True:
    # Bed file directories
    mod_bed_file_path = os.path.join(ref_dir, mod_base + '.bed')
    can_bed_file_path = os.path.join(ref_dir, can_base + '.bed')

    # Read the .bed file into a Pandas DataFrame
    mod_bed_df = pd.read_csv(mod_bed_file_path, delimiter=r'\s+', header=None)
    can_bed_df = pd.read_csv(can_bed_file_path, delimiter=r'\s+', header=None)

    # Name the columns for better readability
    mod_bed_df.columns = ['Alignment', 'Start', 'End', 'Name', 'Score', 'Strand']
    can_bed_df.columns = ['Alignment', 'Start', 'End', 'Name', 'Score', 'Strand']

    # Filter the DataFrame to only keep rows where 'Alignment' equals the specified alignments
    mod_bed_df = mod_bed_df.query(f'Alignment == "{mod_alignment}"')
    can_bed_df = can_bed_df.query(f'Alignment == "{can_alignment}"')
    #print(mod_bed_df)
    #print(can_bed_df)
    mod_bed_df.to_csv(mod_bed_file_path, sep='\t', header=False, index=False)
    can_bed_df.to_csv(can_bed_file_path, sep='\t', header=False, index=False)


if os.stat(os.path.join(ref_dir,mod_base+'.bed')).st_size == 0 or os.stat(os.path.join(ref_dir,can_base+'.bed')).st_size == 0: 
    print('Xemora [ERROR] - Empty bed file generated. Check that XNA bases were present in sequence of input fasta file, and correct XNA base specified. Also check if bed_filtering is set')
    sys.exit()
################################################################################################################################


#      '+os.path.join(mod_bam_dir,os.path.basename(mod_bam_dir))+'.bam'+' \



################################################################################################################################
#Step 5: Generate Chunks. 

#Default name for all bam path directory 
filtered_bam_path = 'bam/bam_filtered.bam'

if regenerate_chunks == True:
    print('Xemora [STATUS] - Generating chunks for modified base training.')
    cmd = 'remora \
      dataset prepare \
      '+os.path.join(mod_pod_dir,os.path.basename(xna_raw_dir))+'.pod5'+' \
       '+os.path.join(mod_dir, filtered_bam_path) + ' \
      --output-remora-training-file '+os.path.join(chunk_dir,'mod_chunks.npz')+' \
      --focus-reference-positions '+os.path.join(ref_dir,mod_base)+'.bed'+' \
      --mod-base '+mod_base+' '+mod_base+' \
      --motif '+can_base+' 0 \
      --kmer-context-bases '+kmer_context+' \
      --refine-kmer-level-table '+kmer_table_path+' \
      --refine-rough-rescale '+' \
      --chunk-context '+chunk_context
    os.system(cmd)

#      --focus-reference-positions '+os.path.splitext(os.path.join(ref_dir,'x'+os.path.basename(xna_ref_fasta)))[0]+'.bed'+' \
#      '+os.path.join(can_bam_dir,os.path.basename(can_bam_dir))+'.bam'+' \
if regenerate_chunks == True: 
    print('Xemora [STATUS] - Generating chunks for canonical base training.')
    cmd = 'remora \
      dataset prepare \
      '+os.path.join(can_pod_dir,os.path.basename(dna_raw_dir))+'.pod5'+' \
       '+os.path.join(can_dir, filtered_bam_path) + ' \
      --output-remora-training-file '+os.path.join(chunk_dir,'can_chunks.npz')+' \
      --focus-reference-positions '+os.path.join(ref_dir,can_base)+'.bed'+' \
      --mod-base-control \
      --motif '+can_base+' 0 \
      --kmer-context-bases '+kmer_context+' \
      --refine-kmer-level-table '+kmer_table_path+' \
      --refine-rough-rescale '+' \
      --chunk-context '+chunk_context
    os.system(cmd)
################################################################################################################################



if remerge_chunks == True: 
    print('Xemora [STATUS] - Merging chunks for training.')
    cmd = 'remora \
      dataset merge \
      --balance \
      --input-dataset '+os.path.join(chunk_dir,'mod_chunks.npz')+' '+chunk_num+'_000 \
      --input-dataset '+os.path.join(chunk_dir,'can_chunks.npz')+' '+chunk_num+'_000 \
      --output-dataset '+os.path.join(chunk_dir,'training_chunks.npz')
    os.system(cmd)

if gen_model == True:
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
      --batch-size 100 '# + '\
      #--ext-val ' + ext_dir

    os.system(cmd)



    print('Xemora [Status] - Complete. Saving model to '+model_dir)


