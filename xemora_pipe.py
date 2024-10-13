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

working_dir = '/home/xenolab/DataAnalysis/Kaplan/training/10.4.1/BSn/241002_BA-ST_REV/ST-REV' 
xna_fast5_dir = '/home/xenolab/DataAnalysis/Kaplan/raw/240627_NTC_Phusion_xr_Train/20240627_1621_MN37138_AUD804_5ac1717b/750_pod5_train'
xna_ref_fasta = '/home/xenolab/DataAnalysis/Kaplan/raw/240627_NTC_Phusion_xr_Train/reference/GAC_1F_7R.fasta'
dna_fast5_dir = '/home/xenolab/DataAnalysis/Kaplan/raw/240104_BSn_90mer_xr_train/20240104_1448_MN37138_ARV509_33c529d5/300_fast5_training'
dna_ref_fasta = '/home/xenolab/DataAnalysis/Kaplan/raw/240104_BSn_90mer_xr_train/BS_90mer.fasta'

############################################################
#Basecall paths

bc_working_dir = '/home/xenolab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/241011_B18_B19_part1_bc/ST-Basecall'
bc_fast5_dir = '/home/xenolab/DataAnalysis/Kaplan/raw/241011_B18_B19_Pol_Screen_1/20241011_1215_MN41475_AVY513_81500ab7/pod5'
bc_xna_ref_fasta = '/home/xenolab/DataAnalysis/Kaplan/raw/241011_B18_B19_Pol_Screen_1/reference/B18_B19_ref.fasta'
bc_model_file = '/home/xenolab/github/kaplanna/xemora/models/240930_NTC_Models/GBC-ST-model_best.pt'
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
