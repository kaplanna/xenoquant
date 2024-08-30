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

working_dir = '/home/marchandlab/DataAnalysis/Kaplan/training/10.4.1/BSn/240627_TSC_Extension_rerun/ST-Train_Q7_val_prop_test_JS'
xna_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/240701_TSC_Extension_xr_Train_rerun/20240701_1546_MN37138_AUD665_f194304a/300_pod5'
xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/240627_TSC_Extension_xr_Train/reference/TSC_Extension_Ref.fasta'
dna_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/240627_NTC_Phusion_xr_Train/20240627_1621_MN37138_AUD804_5ac1717b/800_pod5'
dna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/240627_NTC_Phusion_xr_Train/reference/NTC_Phusion_Ref.fasta'
############################################################
#Basecall paths

bc_working_dir = '/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/240819_B13_Full/BA-Basecall'
bc_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/240816_B13_pH_Rerun_2/20240816_1605_MN41475_AUB839_d3996fad/pod5'
bc_xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/240816_B13_pH_Rerun_2/reference/B13.fasta'
bc_model_file = '/home/marchandlab/DataAnalysis/Kaplan/training/240312_BS_Models_Xenoffice/BA_Train_Q7/model/model_best.pt'
#bc_model_file = working_dir+'/model/model_best.pt'


############################################################

############################################################
train_model = False
basecall_reads = True
output_basecall_results = True
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
    
