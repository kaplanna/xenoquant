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

working_dir = '/home/marchandlab/DataAnalysis/Kaplan/training/10.4.1/BSn/241002_BA-ST_REV/BA-REV' 
xna_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/240927_GGC_xr_90mer/20240927_1237_MN41475_AVW803_5af5c123/2_pod5_train'
xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/240927_GGC_xr_90mer/reference/GGC_12F_20R.fasta'
dna_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/240104_BSn_90mer_xr_train/20240104_1448_MN37138_ARV509_33c529d5/300_fast5_train'
dna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/240104_BSn_90mer_xr_train/BS_90mer.fasta'

############################################################
#Basecall paths

bc_working_dir = '/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/241002_Model_Testing/GBC_Context/XNA_Dataset_90mers/BA-Basecall'
bc_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/240104_BSn_90mer_xr_train/20240104_1448_MN37138_ARV509_33c529d5/200_fast5_test'
bc_xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/240104_BSn_90mer_xr_train/BS_90mer.fasta'
bc_model_file = '/home/marchandlab/github/kaplanna/xemora/models/240930_NTC_Models/GBC-BA-model_best.pt'
#bc_model_file = working_dir+'/model/model_best.pt'


############################################################

############################################################
train_model = False
basecall_reads = True
output_alignment_results = True
cutadapt_demux = False
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
