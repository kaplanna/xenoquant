########################################################################
########################################################################
"""
xr_params.py 

Title: Unpublished work

By: H. Kawabe, N. Kaplan, J. A. Marchand

Updated: 3/2/23
"""
########################################################################
########################################################################

import numpy as np



#Standard basepairs written in 'purine pyrimidine' order
standard_base_pairs = ['AT','GC']

#Convert this to set
standard_bases = np.concatenate(list(list(i) for i in standard_base_pairs))

#Alternative basepairs written in 'purine pyrimidine' order
xna_base_pairs = ['BS','PZ','JV','XK']

#Specify canonical base substitution desired for xFASTA generation here
confounding_pairs =  ['BA','SA','PG','ZC','JC','VG','XA','KG'] 

#If XNAs are given different standard base substitutions, set them up as seperate (e.g, ['P','Z'])
xna_segmentation_model_sets = ['BS','PZ','JV','XK', 'QW','ER']

#Possible XNA bases
xna_bases = np.concatenate(list(list(i) for i in xna_base_pairs))


######################XFASTA GENERATION######################

#Fasta2x - write sequences to xfasta even if they do no contain XNAs. Default = False 
write_no_xna_seq = True

#Fasta2x - Write gaps in place of XNAs in fasta reference file for null testing
write_gaps = False

############################################################
##Analysis instructions 

#Re-basecall pod5 file. Required if new reference files are being used. 
basecall_pod = False

#Re-generate BAM files for reference-based basecalling.
regenerate_bam = False

#Re-generate training or basecalling chunks.
regenerate_chunks = True
############################################################


############################################################
##Model Training and Basecalling Parameters

#Modified base in Fasta sequence you wish to train model or use model to basecall
mod_base = 'Z'

#Most similar substituted canonical base you will be comparing against 
can_base = 'C'

#Extent of Kmer content (-,+) to store for model training
kmer_context ='4 4' 

#Extent of chunk context (centered around modified base) 
chunk_context = '50 50' 


############################################################
#Guppy Base caller configuration

#Path to guppy basecaller
basecaller_path ='~/ont-guppy/bin/guppy_basecaller' 

#GPU enabled 
device_type = 'cuda:all' 

#Config file 
guppy_config_file = 'dna_r9.4.1_450bps_hac.cfg' 


