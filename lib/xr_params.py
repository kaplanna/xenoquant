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
standard_base_pairs = ['AT','GC', 'NN']

#Convert this to set
standard_bases = np.concatenate(list(list(i) for i in standard_base_pairs))

#Alternative basepairs written in 'purine pyrimidine' order
xna_base_pairs = ['BS','PZ','JV','XK']

#Specify canonical base substitution desired for xFASTA generation here
confounding_pairs =  ['BA','SA','PG','ZC','JC','VG','XA','KG'] 

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

#Re-basecall pod5 file. Required if new reference files are being used. 
basecall_pod = True

#Dual barcode basecall
barcode_basecall = True


#Re-generate BAM files for reference-based basecalling.
regenerate_bam = True

#Converting BAM files for data correction 
bam_convert = False

#Merge fail bam ?
merge_fail = False

#convert bam files to fasta for troubleshooting- NOTE: use xemora-re env
bam_to_fasta = False

#Filtering bed files by reference sequence - only use if training on mixed data sets
bed_filtering = True
mod_alignment = "P_Xemora_Train+XPOS[P:130]"
can_alignment = "G_Xemora_Train+XPOS[P:130]"

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
mod_base = 'Z'

#Most similar substituted canonical base you will be comparing against 
can_base = 'C'

#Extent of Kmer content (-,+) to store for model training
kmer_context ='4 4' 

#Extent of chunk context (centered around modified base) 
chunk_context = '50 50' 

#Proportion of reads to use for validation 
val_proportion = '0.2'

#Number of chunks for training (in thousands: e.g.: '200' = 200,000 chunks) 
chunk_num = '500000'

############################################################
# NanoPlot QC Analysis
NanoPlot_Training = False
NanoPlot_Basecall = False


############################################################
#Guppy Base caller configuration

#Path to guppy basecaller
basecaller_path ='~/ont-guppy/bin/guppy_basecaller' 
guppy_barcoder_path ='~/ont-guppy/bin/guppy_barcoder'

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

        





