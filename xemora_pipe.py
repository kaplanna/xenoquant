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

working_dir = '/home/xenolab/DataAnalysis/Kaplan/training/10.4.1/BSn/240930_NTC_Training_Set/240930_GBT_Train/ST-Train' 
xna_fast5_dir = '/home/xenolab/DataAnalysis/Kaplan/raw/240417_ACGBTCG_xr_train/20240417_1649_MN37138_ASP586_df7c41f0/200_pod5_training'
xna_ref_fasta = '/home/xenolab/DataAnalysis/Kaplan/raw/240417_ACGBTCG_xr_train/reference/ACGBTCG_ref.fasta'
dna_fast5_dir = '/home/xenolab/DataAnalysis/Kaplan/raw/240627_NTC_Phusion_xr_Train/20240627_1621_MN37138_AUD804_5ac1717b/750_pod5_train'
dna_ref_fasta = '/home/xenolab/DataAnalysis/Kaplan/raw/240627_NTC_Phusion_xr_Train/reference/GAT_5F_10R.fasta'

############################################################
#Basecall paths

bc_working_dir = '/Users/nickkaplan/xenobiocode/DataAnalysis/basecall/240830_B16_B17_Basecall/BA-Basecall'
bc_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/240816_B13_pH_Rerun_2/20240816_1605_MN41475_AUB839_d3996fad/pod5'
bc_xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/240816_B13_pH_Rerun_2/reference/B13.fasta'
bc_model_file = '/home/marchandlab/DataAnalysis/Kaplan/training/240312_BS_Models_Xenoffice/BA_Train_Q7/model/model_best.pt'
#bc_model_file = working_dir+'/model/model_best.pt'


############################################################

############################################################
train_model = True
basecall_reads = False
output_basecall_results = False
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
if output_basecall_results==True: 
    script_path = './lib/xr_results.py'
    cmd = f'python {script_path} {bc_working_dir}'
    os.system(cmd)
    

#demux
if cutadapt_demux==True: 
    demux_path = './demux/xr_demux_cutadapt.py'
    cmd = f'python {demux_path} {bc_working_dir}'
    os.system(cmd)    
