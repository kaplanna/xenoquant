########################################################################
########################################################################
"""
xr_train.py 

Title: Unpublished work

By: H. Kawabe, N. Kaplan, J. A. Marchand

Updated: 3/2/23
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
xna_fast5_dir = os.path.expanduser(sys.argv[2])
xna_ref_fasta = os.path.expanduser(sys.argv[3])
dna_fast5_dir = os.path.expanduser(sys.argv[4])
dna_ref_fasta = os.path.expanduser(sys.argv[5])

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
    print('Xemora [ERROR] - Reference fasta file not file. Please check file exist or file path.')
    sys.exit()


#Step 1: Generate pod5 files for modified base
if os.path.isfile(os.path.join(mod_pod_dir,os.path.basename(xna_fast5_dir))+'.pod5')==False: 
    cod5_to_fast5(get_fast5_subdir(xna_fast5_dir), os.path.join(mod_pod_dir,os.path.basename(xna_fast5_dir))+'.pod5')
else: 
    print('Xemora [STATUS] - POD5 file for modified base found. Skipping POD5 coversion')


if os.path.isfile(os.path.join(can_pod_dir,os.path.basename(dna_fast5_dir))+'.pod5')==False: 
    cod5_to_fast5(get_fast5_subdir(dna_fast5_dir), os.path.join(can_pod_dir,os.path.basename(dna_fast5_dir))+'.pod5')
else: 
    print('Xemora [STATUS] - POD5 file for canonical base found. Skipping POD5 coversion')


#Step 2: #Basecall pod5 files 
if basecall_pod ==True: 
    cmd=os.path.expanduser(basecaller_path)+' -i '+mod_pod_dir+' -s '+mod_fastq_dir+' -c '+guppy_config_file+' -x auto --bam_out --index --moves_out -a '+os.path.join(ref_dir,'x'+os.path.basename(xna_ref_fasta))
    os.system(cmd)
else: 
    print('Xemora [STATUS] - Skipping POD5 basecalling for modified bases.')


if basecall_pod ==True: 
    cmd=os.path.expanduser(basecaller_path)+' -i '+can_pod_dir+' -s '+can_fastq_dir+' -c '+guppy_config_file+' -x auto --bam_out --index --moves_out -a '+os.path.join(ref_dir,'x'+os.path.basename(dna_ref_fasta))
    os.system(cmd)
else: 
    print('Xemora [STATUS] - Skipping POD5 basecalling for canonical bases.')



#Step 3: Merge Bam files 
if os.path.isfile(os.path.join(mod_bam_dir,os.path.basename(mod_bam_dir))+'.bam') == False or regenerate_bam == True: 
    cmd = 'samtools merge '+os.path.join(mod_bam_dir,os.path.basename(mod_bam_dir))+'.bam'+' '+os.path.join(mod_fastq_dir,'pass/*.bam -f')
    print('Xemora [STATUS] - Merging modified BAM files.')
    os.system(cmd)

if os.path.isfile(os.path.join(can_bam_dir,os.path.basename(can_bam_dir))+'.bam') == False or regenerate_bam == True: 
    cmd = 'samtools merge '+os.path.join(can_bam_dir,os.path.basename(can_bam_dir))+'.bam'+' '+os.path.join(can_fastq_dir,'pass/*.bam -f')
    print('Xemora [STATUS] - Merging canonical BAM files.')
    os.system(cmd)




#Step 4: Bed file generation 
if os.stat(os.path.join(ref_dir,'x'+os.path.basename(xna_ref_fasta))).st_size == 0: 
    print('Xemora [ERROR] - Empty xfasta file generated. Check that XNA bases were present in sequence of input fasta file.')
    sys.exit()

print('Xemora [STATUS] - Generating bed file for modified base.')
cmd = 'python lib/xr_xfasta2bed.py '+os.path.join(ref_dir,'x'+os.path.basename(xna_ref_fasta))+' '+os.path.join(ref_dir,mod_base+'.bed ' +mod_base+' '+mod_base)
os.system(cmd)

print('Xemora [STATUS] - Generating bed file for canonical base.')
cmd = 'python lib/xr_xfasta2bed.py '+os.path.join(ref_dir,'x'+os.path.basename(xna_ref_fasta))+' '+os.path.join(ref_dir,can_base+'.bed '+mod_base+' '+can_base)
os.system(cmd)



if os.stat(os.path.join(ref_dir,mod_base+'.bed')).st_size == 0 or os.stat(os.path.join(ref_dir,can_base+'.bed')).st_size == 0: 
    print('Xemora [ERROR] - Empty bed file generated. Check that XNA bases were present in sequence of input fasta file, and correct XNA base specified.')
    sys.exit()



#Step 5: Generate Chunks. 
if regenerate_chunks == True: 
    print('Xemora [STATUS] - Generating chunks for modified base training.')
    cmd = 'remora \
      dataset prepare \
      '+os.path.join(mod_pod_dir,os.path.basename(xna_fast5_dir))+'.pod5'+' \
      '+os.path.join(mod_bam_dir,os.path.basename(mod_bam_dir))+'.bam'+' \
      --output-remora-training-file '+os.path.join(chunk_dir,'mod_chunks.npz')+' \
      --focus-reference-positions '+os.path.join(ref_dir,mod_base)+'.bed'+' \
      --mod-base '+mod_base+' '+mod_base+' \
      --motif '+can_base+' 0 \
      --kmer-context-bases '+kmer_context+' \
      --chunk-context '+chunk_context
    os.system(cmd)

#      --focus-reference-positions '+os.path.splitext(os.path.join(ref_dir,'x'+os.path.basename(xna_ref_fasta)))[0]+'.bed'+' \

if regenerate_chunks == True: 
    print('Xemora [STATUS] - Generating chunks for canonical base training.')
    cmd = 'remora \
      dataset prepare \
      '+os.path.join(can_pod_dir,os.path.basename(dna_fast5_dir))+'.pod5'+' \
      '+os.path.join(can_bam_dir,os.path.basename(can_bam_dir))+'.bam'+' \
      --output-remora-training-file '+os.path.join(chunk_dir,'can_chunks.npz')+' \
      --focus-reference-positions '+os.path.join(ref_dir,can_base)+'.bed'+' \
      --mod-base-control \
      --motif '+can_base+' 0 \
      --kmer-context-bases '+kmer_context+' \
      --chunk-context '+chunk_context
    os.system(cmd)

if regenerate_chunks == True: 
    print('Xemora [STATUS] - Merging chunks for training.')
    cmd = 'remora \
      dataset merge \
      --input-dataset '+os.path.join(chunk_dir,'mod_chunks.npz')+' 40_000 \
      --input-dataset '+os.path.join(chunk_dir,'can_chunks.npz')+' 40_000 \
      --output-dataset '+os.path.join(chunk_dir,'training_chunks.npz')
    os.system(cmd)


print('Xemora [STATUS] - Training model.')
cmd = 'remora \
  model train \
  '+os.path.join(chunk_dir,'training_chunks.npz')+' \
  --model models/Conv_w_ref.py \
  --device 0 \
  --output-path '+model_dir+' \
  --overwrite \
  --kmer-context-bases '+kmer_context+' \
  --chunk-context '+chunk_context+' \
  --batch-size 100'
os.system(cmd)



print('Xemora [Status] - Complete. Saving model to '+model_dir)


