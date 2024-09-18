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

#Specify canonical base substitution desired for xFASTA generation here #changed to ST by NK
confounding_pairs =  ['BA','ST','PG','ZC','JC','VG','XA','KG'] 

#If XNAs are given different standard base substitutions, set them up as seperate (e.g, ['P','Z'])
xna_segmentation_model_sets = ['B','S','PZ','JV','X', 'K', 'QW','ER']

#Possible XNA bases
xna_bases = np.concatenate(list(list(i) for i in xna_base_pairs))


######################XFASTA GENERATION######################

#Fasta2x - write sequences to xfasta even if they do no contain XNAs. Default = False 
write_no_xna_seq = False

#Fasta2x - Write gaps in place of XNAs in fasta reference file for null testing
write_gaps = False

############################################################
##Analysis instructions 

#Re-basecall pod5 file. Required if new reference files are being used. CHANGE
basecall_pod = True

#Dual barcode basecall
barcode_basecall = False


#Re-generate BAM files for reference-based basecalling. CHANGE
regenerate_bam = True

#Converting BAM files for data correction 
bam_convert = False

#Merge fail bam ?
merge_fail = False

#convert bam files to fasta for troubleshooting- NOTE: use xemora-re env
bam_to_fasta = False

#Filtering bed files by reference sequence - only use if training on mixed data sets
bed_filtering = False
mod_alignment = "TSC_Extended+XPOS[B:105]"
can_alignment = "TTC_4F_9R+XPOS[B:75]"

#Data extraction, filtering, and heptamer correction 
data_fix = True

#SAM file sequence corrections 
sam_corrections = True

#Convering SAM files for training 
sam_convert = True

#Re-generate training or basecalling chunks.
regenerate_chunks = True

#Merge chunks again for training data. 
remerge_chunks = True

#Build model using Remora 
gen_model = True
############################################################


############################################################
##Model Training and Basecalling Parameters

#kmer table 
#kmer_table_path = 'models/remora/4mer_9.4.1.csv'
kmer_table_path = 'models/remora/9mer_10.4.1.csv'

#ml model (ConvLSTM_w_ref.py or Conv_w_ref.py')
ml_model_path = 'models/ConvLSTM_w_ref.py'


#Modified base in Fasta sequence you wish to train model or use model to basecall
mod_base = 'B'

#Most similar substituted canonical base you will be comparing against 
can_base = 'A'

#Extent of Kmer content (-,+) to store for model training
kmer_context ='4 4' 

#Extent of chunk context (centered around modified base) 
chunk_context = '50 50' 

#Proportion of reads to use for validation 
val_proportion = '0.2'

#Number of chunks for training (in thousands: e.g.: '200' = 200,000 chunks) 
chunk_num = '500000'




############################################################
# New parameters for xr_train_methods.py from Jayson
overwrite_pod = True
dorado_path = '~/dorado-0.7.2-linux-x64/bin/dorado'
dorado_model = '~/dorado-0.7.2-linux-x64/bin/dna_r10.4.1_e8.2_400bps_hac@v5.0.0'
min_qscore = 7
#Range of chunk context to use (in bp) for modified base training (default +/- 0) 
mod_chunk_range = 0
can_chunk_range = 0

#Shift the mod chunk range position by a fixed amount (default = 0) 
mod_chunk_shift = 0
can_chunk_shift = 0

#Balance training chunks. May be set to false for testing, otherwise set to true. 
balance_chunks = True

max_mod_reads = 0
max_can_reads = 0

filter_mod_readIDs = ''
filter_can_readIDs = ''#'/home/xenolab/DataAnalysis/Kaplan/basecall/240627_NTC_Phusion_xr_Train_Basecall/cutadapt_demux/NB02_FWD_NB08_REV_read_ids.txt'

############################################################
# NanoPlot QC Analysis
NanoPlot_Training = False
NanoPlot_Basecall = False


############################################################
#Guppy Base caller configuration

#Path to guppy basecaller
basecaller_path ='~/ont-guppy/bin/guppy_basecaller' 
guppy_barcoder_path ='~/ont-guppy/bin/guppy_barcoder'
guppy_aligner_path = '~/ont-guppy/bin/guppy_aligner' 

#GPU enabled 
device_type = 'cuda:all' 

#Guppy q-score threshold for pass/fail 
training_mod_min_qscore = 7
training_can_min_qscore = 7
basecall_min_qscore = 7

#Config file 
#guppy_config_file = 'dna_r9.4.1_450bps_hac.cfg'
guppy_config_file = 'dna_r10.4.1_e8.2_400bps_hac.cfg'
#guppy_config_file = 'dna_r10.4.1_e8.2_260bps_hac.cfg'
#guppy_config_file = 'dna_r10.4.1_e8.2_260bps_sup.cfg'

#barcoding
#for dual
#barcode_config = 'configuration_dual.cfg'
#barcode_kit = 'EXP-DUAL00'

       




