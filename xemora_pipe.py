########################################################################
########################################################################
"""
xemora_pipe.py 


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

working_dir = '/home/marchandlab/DataAnalysis/Kaplan/training/10.4.1/BSn/241127_GBG_Training/BA-Train' 
xna_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/241114_GBG-CSC_90mer_xr_train/20241114_1355_MN41475_AWO394_fc4c9c6b/pod5_train'
xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/241114_GBG-CSC_90mer_xr_train/reference/GBG_90mer_ref.fasta'
dna_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/240627_NTC_Phusion_xr_Train/20240627_1621_MN37138_AUD804_5ac1717b/pod5'
dna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/240627_NTC_Phusion_xr_Train/reference/GAG_Ref_NTC.ref'

############################################################
#Basecall paths

bc_working_dir = '/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/XPCR/241126_B26-B27_Basecall/GBG-ST-Basecall'
bc_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/241125_B26_B27_Seq_Context/20241125_1645_MN41475_AWF292_cfbb6444/pod5'
bc_xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/241125_B26_B27_Seq_Context/reference/B26_B27_ref.fasta'
bc_model_file = '/home/marchandlab/github/kaplanna/xemora/models/240930_NTC_Models/GBG-ST-model_best.pt'
#bc_model_file = working_dir+'/model/model_best.pt'


############################################################
############################################################
train_model = False
basecall_reads = True
output_alignment_results = True
cutadapt_demux = True
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
if output_alignment_results==True: 
    results_path = './lib/xr_results.py'
    cmd = f'python {results_path} {bc_working_dir}'
    os.system(cmd)
    
#demux
if cutadapt_demux==True: 
    demux_path = './lib/xr_demux.py'
    cmd = f'python {demux_path} {bc_working_dir}'
    os.system(cmd)    
