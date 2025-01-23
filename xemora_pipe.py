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

working_dir = '/home/marchandlab/DataAnalysis/Kaplan/training/10.4.1/BSn/' 
xna_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/BS_xr/240925_GBA-TSC_anneal_xr_train/20240925_1500_MN37138_AUC702_eb4d5541/300_pod5_train'
xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/BS_xr/240925_GBA-TSC_anneal_xr_train/reference/ref_GBA.fasta'
dna_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/BS_xr/240627_NTC_Phusion_xr_Train/20240627_1621_MN37138_AUD804_5ac1717b/pod5'
dna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/BS_xr/240627_NTC_Phusion_xr_Train/reference/ref_Phusion_GBA.fasta'

############################################################
#Basecall paths

bc_working_dir = '/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/XPCR/250111_B32_B33_Nucleotide_Basecall/ST_Basecall'
bc_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/xPCR/250110_B32_B33_Nucleotide_Expmt_Rerun/20250110_1620_MN41475_AWI420_f963a31c/pod5'
bc_xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/xPCR/250110_B32_B33_Nucleotide_Expmt_Rerun/reference/ref_B32_B33.fasta'
bc_model_file = '/home/marchandlab/github/kaplanna/xemora/models/240930_NTC_Models/GBC-ST-model_best.pt'
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
