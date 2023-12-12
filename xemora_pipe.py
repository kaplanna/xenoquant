########################################################################
########################################################################
"""
xemora_pipe.py 

Title: Unpublished work

By: H. Kawabe, N. Kaplan, J. A. Marchand

Updated: 8/20/23
"""
########################################################################
########################################################################


import os
import numpy as np
import glob
import sys
from pathlib import Path
from lib.xr_tools import *
from lib.xr_params import *


############################################################
#Training paths
'''
working_dir = '/home/marchandlab/DataAnalysis/Kaplan/training/10.4.1/231205_train_230723_Old_PZ/ZC_Train_Q9_w_Fail'
xna_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/fast5/10.4.1/230712_PZ_Xmra_Train_10_4_1/50fast5'
xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/fast5/10.4.1/230712_PZ_Xmra_Train_10_4_1/20230712_1805_MN37138_AOX743_3c1cdb19/reference/P_xr_train_full_adapters.fa'
dna_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/fast5/10.4.1/230712_GC_Xmra_Train_10_4_1/50fast5'
dna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/fast5/10.4.1/230712_GC_Xmra_Train_10_4_1/reference/GC_xem_adapters.fa'
'''

#working_dir = '/home/marchandlab/DataAnalysis/Kawabe/231210_BSn_IDT/xemora/isoGmod_isoCcan'
#xna_fast5_dir = '/home/marchandlab/DataAnalysis/Kawabe/231210_BSn_IDT/20231207_1121_MN41475_AQQ090_fa649842/small_fast5'
#xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kawabe/231210_BSn_IDT/reference/isoG.fasta'

#dna_fast5_dir = '/home/marchandlab/DataAnalysis/Kawabe/231210_BSn_IDT/20231207_1121_MN41475_AQQ090_fa649842/small_fast5'
#dna_ref_fasta = '/home/marchandlab/DataAnalysis/Kawabe/231210_BSn_IDT/reference/isoC.fasta'

working_dir = '/home/marchandlab/DataAnalysis/Kaplan/training/10.4.1/231211_GAP_Train/ZC_Train'
xna_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/fast5/10.4.1/230928_GAP_Train/GAP_Train_fast5_0-49'
xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/fast5/10.4.1/230928_GAP_Train/reference/GAP_allseq_wadapt.fa'
dna_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/fast5/10.4.1/230928_GAP_Train/GAP_Train_fast5_0-49'
dna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/fast5/10.4.1/230928_GAP_Train/reference/GAP_allseq_wadapt.fa'

############################################################
#Basecall paths

bc_working_dir = '/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/231207_barcoding_testing_GAP/ZC_Basecall'
bc_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/fast5/10.4.1/230928_GAP_Train/GAP_basecall_92fast5'
bc_xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/fast5/10.4.1/230928_GAP_Train/reference/GAP_allseq_wadapt.fa'
#bc_model_file = '/home/marchandlab/DataAnalysis/Kaplan/training/10.4.1/230928_GAP_Train/ZC_Train/model/model_best.pt'
bc_model_file = working_dir+'/model/model_best.pt'


############################################################

############################################################
train_model = True
basecall_reads = False
output_basecall_results = False
############################################################

#conda activate xemora-re


        
#Train dataset with xemora train
if train_model ==True: 
    cmd = 'python xemora.py train -w '+working_dir+' -f '+xna_fast5_dir+' '+dna_fast5_dir+' -r '+xna_ref_fasta+' '+dna_ref_fasta
    os.system(cmd)


#Basecall fast5 directory 
if basecall_reads==True: 
    cmd = 'python xemora.py basecall -w '+bc_working_dir+' -f '+bc_fast5_dir+' -r '+bc_xna_ref_fasta+' -m '+bc_model_file 
    os.system(cmd)


#output results
if output_basecall_results==True: 
    script_path = './lib/xr_results.py'
    cmd = f'python {script_path} {bc_working_dir}'
    os.system(cmd)
    
