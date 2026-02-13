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

working_dir = '/home/marchandlab/DataAnalysis/Kaplan/training/2509_Signal_Plots/250922_CNT_plots/Px-N'
xna_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/DsPx/241209_C_Ds_T_Px_diol/20241209_1730_MN37138_AWF971_532b88f2/50_pod5_train'
xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/DsPx/241209_C_Ds_T_Px_diol/reference/CDsT_60mer.fasta'
dna_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/DsPx/250714_CNT_xr_train_JS_767/20250714_1733_MN41475_AYG523_b8fd4238/pod5_train'
dna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/DsPx/241209_C_Ds_T_Px_diol/reference/CDsT_60mer.fasta'

############################################################
#Basecall paths

bc_working_dir = '/home/marchandlab/DataAnalysis/Kaplan/basecall/8letter/251126_H7-8_Basecall/PG-Basecall_100-100'
bc_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/xPCR/251125_H7-8_6LT-8TP/20251125_1800_MN37138_AYY072_d282dce0/pod5'
bc_xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/xPCR/251125_H7-8_6LT-8TP/reference/REF_H7-8_P_only.fasta'
barcode_pair_csv = '/home/marchandlab/DataAnalysis/Kaplan/raw/xPCR/251125_H7-8_6LT-8TP/reference/DEMUX_H7-8.csv'
bc_model_file = '/home/marchandlab/DataAnalysis/Kaplan/training/PZ/251020_PG_training/PG-Train_100_100/model/model_best.pt'

#bc_model_file = working_dir+'/model/model_best.pt'


############################################################
############################################################
train_model = False
basecall_reads = True
output_alignment_results = True
cutadapt_demux = True
raw_basecall_analysis = False
##########################################################
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
    cmd = f'python {demux_path} {bc_working_dir} {barcode_pair_csv}'
    os.system(cmd)    
    
#raw basecall analysis
if raw_basecall_analysis==True: 
    raw_bc_path = './lib/xr_raw_basecall_analysis.py'
    cmd = f'python {raw_bc_path} {bc_working_dir} {FILTER_BY_CLASS_0}'
    os.system(cmd)    
