########################################################################
########################################################################
"""
xr_basecall.py 

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
working_dir = os.path.expanduser(sys.argv[1])
fast5_dir = os.path.expanduser(sys.argv[2])
xna_ref_fasta = os.path.expanduser(sys.argv[3])
model_file = os.path.expanduser(sys.argv[4])
############################################################




#Generate directories
print('Xemora [Status] - Initializing Xemora basecalling.')


working_dir = check_make_dir(working_dir)
ref_dir = check_make_dir(os.path.join(working_dir,'references'))
chunk_dir = check_make_dir(os.path.join(working_dir,'chunks'))
mod_dir = check_make_dir(os.path.join(working_dir,'preprocess'))
mod_pod_dir = check_make_dir(os.path.join(mod_dir,'pod5'))
mod_fastq_dir = check_make_dir(os.path.join(mod_dir,'fastq'))
mod_bam_dir = check_make_dir(os.path.join(mod_dir,'bam'))



#Step 0: FASTA to xFASTA conversion
if os.path.isfile(os.path.expanduser(xna_ref_fasta)): 
    cmd = 'python lib/xr_fasta2x_rc.py '+os.path.expanduser(xna_ref_fasta)+' '+os.path.join(ref_dir,'x'+os.path.basename(xna_ref_fasta))
    os.system(cmd)
else: 
    print('Xemora  [ERROR] - Reference fasta xna file not file. Please check file exist or file path.')
    sys.exit()


#Step 1: Generate pod5 files for modified base
if os.path.isfile(os.path.join(mod_pod_dir,os.path.basename(fast5_dir))+'.pod5')==False: 
    cod5_to_fast5(get_fast5_subdir(fast5_dir), os.path.join(mod_pod_dir,os.path.basename(fast5_dir))+'.pod5')
else: 
    print('Xemora  [STATUS] - POD5 file for modified base found. Skipping POD5 coversion')


#Step 2: #Basecall pod5 files 
if basecall_pod ==True: 
    cmd=os.path.expanduser(basecaller_path)+' -i '+mod_pod_dir+' -s '+mod_fastq_dir+' '+guppy_config_file+' -x auto --bam_out --index --moves_out -a '+os.path.join(ref_dir,'x'+os.path.basename(xna_ref_fasta))
    os.system(cmd)
else: 
    print('Xemora  [STATUS] - Skipping POD5 basecalling for modified bases.')


#Step 3: Merge Bam files 
if os.path.isfile(os.path.join(mod_bam_dir,os.path.basename(mod_bam_dir))+'.bam') == False or regenerate_bam == True: 
    cmd = 'samtools merge '+os.path.join(mod_bam_dir,os.path.basename(mod_bam_dir))+'.bam'+' '+os.path.join(mod_fastq_dir,'pass/*.bam -f')
    print('Xemora  [STATUS] - Merging modified BAM files.')
    os.system(cmd)


#Step 4: Bed file generation 
if os.stat(os.path.join(ref_dir,'x'+os.path.basename(xna_ref_fasta))).st_size == 0: 
    print('Xemora  [ERROR] - Empty xfasta file generated. Check that XNA bases were present in sequence of input fasta file.')
    sys.exit()

print('Xemora  [STATUS] - Generating bed file for modified base.')
cmd = 'python lib/xr_xfasta2bed.py '+os.path.join(ref_dir,'x'+os.path.basename(xna_ref_fasta))+' '+os.path.join(ref_dir,mod_base+'.bed ' +mod_base+' '+mod_base)
os.system(cmd)



#Step 5: Generate Chunks. 
if regenerate_chunks == True: 
    print('Xemora  [STATUS] - Generating chunks for modified basecalling.')
    cmd = 'remora \
      dataset prepare \
      '+os.path.join(mod_pod_dir,os.path.basename(fast5_dir))+'.pod5'+' \
      '+os.path.join(mod_bam_dir,os.path.basename(mod_bam_dir))+'.bam'+' \
      --output-remora-training-file '+os.path.join(chunk_dir,'basecall_chunks.npz')+' \
      --focus-reference-positions '+os.path.join(ref_dir,mod_base)+'.bed'+' \
      --mod-base '+mod_base+' '+mod_base+' \
      --motif '+can_base+' 0 \
      --kmer-context-bases '+kmer_context+' \
      --chunk-context '+chunk_context
    os.system(cmd)

try: 
    print('Xemora  [STATUS] - Performing basecalling.')
    cmd = 'remora \
      validate from_remora_dataset \
      '+os.path.join(chunk_dir,'basecall_chunks.npz')+' \
      --model '+os.path.expanduser(model_file)+' \
      --full-results-filename '+os.path.join(working_dir,'per-read_modifications.tsv')+' \
      --out-file '+os.path.join(working_dir,'summary_modifications.tsv')
    os.system(cmd)

    print('Xemora  [STATUS] - Basecalling done.')
    print('Xemora  [STATUS] - Basecalling done. Saving results '+os.path.join(working_dir,'per-read_modifications.tsv'))
    print('Xemora  [STATUS] - Basecalling done. Saving results '+os.path.join(working_dir,'summary_modifications.tsv'))
    print('Xemora  [STATUS] - Exiting')
except:
    print('Xemora  [ERROR] - Failed to initialize basecalling model. Check logs.')

