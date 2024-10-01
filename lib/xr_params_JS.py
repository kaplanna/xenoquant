########################################################################
########################################################################
"""
xr_params.py 

Title: Unpublished work

By: H. Kawabe, N. Kaplan, J. Sumabat, J. A. Marchand

Updated: 9/5/24
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
#xna_base_pairs = ['DP','BS','JV','XK']

#Specify canonical base substitution desired for xFASTA generation here
confounding_pairs =  ['DC','BA','ST','PN','ZN','JC','VG','XA','KG'] 

#If XNAs are given different standard base substitutions, set them up as seperate (e.g, ['P','Z'])
xna_segmentation_model_sets = ['B','S','P', 'Z','JV','X', 'K', 'QW','ER']
#xna_segmentation_model_sets = ['B','S','D','P','JV','X', 'K', 'QW','ER']

#Possible XNA bases
xna_bases = np.concatenate(list(list(i) for i in xna_base_pairs))

############################XFASTA GENERATION###########################
#Fasta2x - write sequences to xfasta even if they do no contain XNAs. Default = False 
write_no_xna_seq = False

#Fasta2x - Write gaps in place of XNAs in fasta reference file for null testing
write_gaps = False

########################################################################
########################SHARED PARAMETERS###############################
##Verbose output for context method operations
verbose_output = False

#Merging and conversion of pod5 files. If true, will overwrite previous pod5 file operations if present. 
overwrite_pod = True

#Re-basecall pod5 file. Required if new reference files are being used. 
basecall_pod = False

#Re-generate BAM files for reference-based basecalling.
regenerate_bam = True

#Re-generate training or basecalling chunks.
regenerate_chunks = True

#Merge chunks again for training data. 
remerge_chunks = False

#Build model using Remora 
gen_model = True

#Display data visuals
display_vis =False

#Run Remora validate and consensus analysis for basecalling
run_analysis = True

#Perform pairwise analysis 
pairwise = False
#################DORADO BASECALLING FILTERS################################

#Max reads to basecall (default = 0 == all reads)
max_reads = 0 

#Q-score filtering threshold for reads during basecalling
min_qscore = 7

#Q-score filtering threshold of basecalled reads (default = 7)
min_qscore_analysis = 7

#File path to set of specific readIDs to basecall; must be new line delimied (default = ''). Will not be used if basecall_by_cluster is set. 
filter_readIDs_bc = ''

########################################################################

############################CONSENSUS RESULTS ANALYSIS#####################
#Minimum amount of unique sequences per reference that had Remora validate done for consensus analysis
min_coverage = 50

#Range of chunk context to use (in bp) for modified base training or basecalling (default +/- 0). For basecalling, If negative will use full read length (e.g., = -1)
mod_chunk_range = 0

#Shift the mod chunk range position by a fixed amount (default = 0) 
mod_chunk_shift = 0


########################################################################
#Basecall reads by vsearch cluster output .txt. Requires reference with matching cluster ids. 
basecall_by_cluster = ''

#Minimum size of cluster to analyzie if cl flag is used for xemora. Note - Cannot be passed as a custom param dictionary. 
min_cluster_size = 100

#Number of clusters to analyze if -cl flag is used for xemora.  Note - Cannot be passed as a custom param dictionary. 
top_n_clusters = 5
########################################################################

##########################MINIMAP 2 Filter##########################
#Basecalling XNA probability filter threshold for LSTM calls (default = 0.5). 
#Unused - For future feature basecall_filter_threshold = 0.5 

#Minimap2 alignment filter for first-round alignment. 60 is top. 0 = no filter. Used in all context training
min_mapq_score = 50

#Minimap2 alignment score filter for first-round alignment. 255 is top. 0?
min_as_score = 130

########################################################################

########################SINGLE CONTEXT TRAINING#########################
#Training only  Modified base in Fasta sequence you wish to train model or use model to basecall
mod_base = 'B'

#Training only  Most similar substituted canonical base you will be comparing against 
can_base = 'A'

##########################ALL CONTEXT TRAINING##########################
#Minimap2 alignment filter for first-round alignment. 60 is top. 0 = no filter. 
min_map_score = 0

########################################################################



#########REMORA CHUNK, TRAINING, & VALIDATION PARAMETERS###############

#kmer table 
#kmer_table_path = 'models/remora/4mer_9.4.1.csv'
#kmer_table_path = '../models/remora/9mer_10.4.1.csv'
kmer_table_path = 'models/remora/9mer_10.4.1.csv'

#ml model (ConvLSTM_w_ref.py or Conv_w_ref.py')
ml_model_path = '../models/ConvLSTM_w_ref.py'
#ml_model_path = 'models/ConvLSTM_w_ref.py'

#Desired Strand
default_strand = '+'

#Number of chunks to include in a batch during model training (default: 1024)
rem_batch = 1024
rem_batch = str(rem_batch)

#Extent of Kmer content (-,+) to store for model training
kmer_context ='4 4' 

#Extent of chunk context (centered around modified base) 
chunk_context = '50 50' 

#Proportion of reads to use for validation 
val_proportion = '0.2'

#Number of chunks for training (in thousands: e.g.: '200' = 200,000 chunks) 
chunk_num = '500000'

########################################################################

############################################################
#Dorado Base caller configuration

#Path to guppy basecaller
dorado_path ='~/dorado-0.8.0-linux-x64/bin/dorado'

#Dorado model file, inputs are either a model path or 'fast', 'hac', or 'sup' for automatic model selection (Default: 'sup')
dorado_model = '/home/marchandlab/dorado-0.8.0-linux-x64/models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0'
#dorado_model = '/home/marchandlab/dorado-0.7.0-linux-x64/models/dna_r10.4.1_e8.2_400bps_hac@v5.0.0'
#dorado_model = '/home/marchandlab/dorado-0.7.0-linux-x64/models/dna_r10.4.1_e8.2_400bps_fast@v5.0.0'
#dorado_model = 'fast'

############################################################
#Guppy Base caller configuration

#Path to guppy basecaller
basecaller_path ='~/ont-guppy/bin/guppy_basecaller' 

#GPU enabled 
device_type = 'cuda:all' 

#Config file 
#guppy_config_file = 'dna_r9.4.1_450bps_hac.cfg'
guppy_config_file = 'dna_r10.4.1_e8.2_400bps_hac.cfg'
#guppy_config_file = 'dna_r10.4.1_e8.2_260bps_hac.cfg'
#guppy_config_file = 'dna_r10.4.1_e8.2_260bps_sup.cfg'
############################################################


##################OLD GUPPY PARAMETER Pipelines##########################
#Merge fail bam (Old parameter, meant for when Guppy basecaller is used)
merge_fail = True

#Data extraction, filtering, and heptamer correction 
data_fix = False

########################################################################
