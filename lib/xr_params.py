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
xna_base_pairs = ['BS','PZ','DX']

#Specify canonical base substitution desired for xFASTA generation here
confounding_pairs =  ['BA','ST','PG','ZC','DA','XT'] 


#If XNAs are given different standard base substitutions, set them up as seperate (e.g, ['P','Z'])
xna_segmentation_model_sets = ['B','S','P','Z','D','X']

#Possible XNA bases
xna_bases = np.concatenate(list(list(i) for i in xna_base_pairs))

require_rc_fasta = True


######################XFASTA GENERATION######################

#Fasta2x - write sequences to xfasta even if they do no contain XNAs. Default = False 
write_no_xna_seq = False

#Fasta2x - Write gaps in place of XNAs in fasta reference file for null testing
write_gaps = False

############################################################
##Analysis instructions 

#Re-basecall pod5 file. Required if new reference files are being used. 
basecall_pod = False

 # Set True to rerun cutadapt demux, False to skip and reuse existing outputs
RERUN_DEMUX = False 

#Re-generate BAM files for reference-based basecalling. DEFAULT: True
regenerate_bam = True

#Converting BAM files for data correction DEFAULT: False
bam_convert = False

#Merge fail bam DEFAULT: False
merge_fail = False

#convert bam files to fasta for troubleshooting- NOTE: use xemora-re env DEFAULT: False
bam_to_fasta = False

#Data extraction, filtering, and heptamer correction 
data_fix = True

#SAM file sequence corrections 
sam_corrections = True

#Convering SAM files for training 
sam_convert = True

#Re-generate training or basecalling chunks.
regenerate_chunks = False

#Generate chunks using basecall anchor (default: False)
bc_anchor = False

#Merge chunks again for training data. 
remerge_chunks = True

#Build model using Remora 
gen_model = False


# ============================
# Alignment filtering
# ============================

# Alignment quality filtering
filter_alignment_fraction = False     # legacy (leave False)
filter_softclip = False

max_total_softclip_frac = 0.30
max_end_softclip_frac   = 0.20
min_aligned_frac        = 0.60


# -------------------------------------------------
# Remora decision threshold override
# -------------------------------------------------

USE_DECISION_THRESHOLD = False     # False = default Remora argmax (0.5)
DECISION_THRESHOLD = 0.90         # Must be > 0.5

############################################################
##Model Training and Basecalling Parameters

#kmer table 
#kmer_table_path = 'models/remora/4mer_9.4.1.csv'

# after (your real layout)
kmer_table_path = "models/remora/9mer_10-4-1.tsv"

#kmer_table_path = 'models/remora/9mer_10-4-1.tsv'
#kmer_table_path = "/home/marchandlab/github/kaplanna/xemora/models/remora/9mer_10-4-1.tsv"



#ml model (ConvLSTM_w_ref.py or Conv_w_ref.py')
ml_model_path = 'models/ConvLSTM_w_ref.py'



#Modified base in Fasta sequence you wish to train model or use model to basecall
mod_base = 'B'
#Most similar substituted canonical base you will be comparing against 
can_base = 'N'

#Extent of Kmer content (-,+) to store for model training
kmer_context ='4 4' #default 4 4 

#Extent of chunk context (centered around modified base) 
chunk_context = '50 50' 

#Proportion of reads to use for validation 
val_proportion = '0.2'

#Number of chunks for training (in thousands: e.g.: '200' = 200,000 chunks) 
chunk_num = '500000' 

############################################################
#Signal Plots

generate_remora_plots = False  # Set to False to skip Remora ref_region plotting

#reference bed flank size for plotting bed file
FLANK = 6

# Always black for canonical
COL_STD = "#000000"

# Default modified color (can override per-project if you want)
COL_MOD = "#e94530" #PZ red
#COL_MOD = "#67BFEB" #BS blue
#COL_MOD = "#ea9f20" #PZ ORANGE

# Highlight fill for XNA (overrideable hex code)
COL_X_HIGHLIGHT = "#9E9E9E"

############################################################
#Raw basecall analysis filter by class 0
FILTER_BY_CLASS_0 = True

############################################################
#CutAdapt Demux parameters
error_rate = 0.2 # default 0.2
min_overlap = 14 # default 14
min_len = 120 # default 120
max_len = 200 # default 200


############################################################
# New parameters for xr_train_methods.py from Jayson
overwrite_pod = False
dorado_path = '~/dorado-0.8.0-linux-x64/bin/dorado'
dorado_model = '~/dorado-0.8.0-linux-x64/models/dna_r10.4.1_e8.2_400bps_hac@v5.0.0'
min_qscore = 5 #default 5
#Range of chunk context to use for modified base training (default +/- 0) 
mod_chunk_range = 0
can_chunk_range = 0

#Shift the mod chunk range position by a fixed amount (default = 0) 
mod_chunk_shift = 0
can_chunk_shift = 0

#Balance training chunks.  May be set to false for testing, otherwise set to true. 
balance_chunks = True

max_mod_reads = 0
max_can_reads = 0
max_bc_reads = 0



#filter_mod_readIDs = '/home/marchandlab/DataAnalysis/Kaplan/raw/8letter/250731_8Letter_PB_GA_Ext/High-conf_training/8L-PG_H6_Highconf_reads_99.txt'
#filter_mod_readIDs = '/home/marchandlab/DataAnalysis/Kaplan/raw/8letter/251210_8L_MixedPZ_train/demux_list/NB25_NB24_read_ids.txt'
filter_mod_readIDs = ''
filter_can_readIDs = '/home/marchandlab/DataAnalysis/Kaplan/training/2509_Signal_Plots/251226_GDsA_AT_Plots/can_set_demux/NB12_FWD_NB20_REV_read_ids.txt'
#filter_can_readIDs = '/home/marchandlab/DataAnalysis/Kaplan/training/BSn/250831_GNT_TNC_training_demux/training_v1_800pod5/demux/NB13_FWD_NB11_REV_read_ids.txt' #for TNC
#filter_readIDs_bc = '/home/marchandlab/DataAnalysis/Kaplan/raw/8letter/250731_8Letter_PB_GA_Ext/251203_8LT_training_demux/NB00_FWD_NB24_REV_read_ids.txt'
filter_readIDs_bc = ''
############################################################


#GPU enabled 
device_type = 'cuda:all' 



#Config file 
#guppy_config_file = 'dna_r9.4.1_450bps_hac.cfg'
guppy_config_file = 'dna_r10.4.1_e8.2_400bps_hac.cfg'
#guppy_config_file = 'dna_r10.4.1_e8.2_260bps_hac.cfg'
#guppy_config_file = 'dna_r10.4.1_e8.2_260bps_sup.cfg'

#barcoding
#for dual
#barcode_config = 'configuration_dual.cfg'
#barcode_kit = 'EXP-DUAL00'

       




