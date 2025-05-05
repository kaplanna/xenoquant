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

working_dir = '/home/marchandlab/DataAnalysis/Kaplan/training/10.4.1/BSn/CBG_xr_train/ST-Train'
xna_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/BS_xr/250408_CBG_xr_train/20250408_1529_MN37138_AYZ259_99a19ebc/400_pod5_train'
xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/BS_xr/250408_CBG_xr_train/reference/ref_CBG.fasta'
dna_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/BS_xr/250408_CAG_xr_train/20250408_1541_MN41475_AYS317_7e48efeb/2_pod5_train'
dna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/BS_xr/250408_CAG_xr_train/20250408_1541_MN41475_AYS317_7e48efeb/reference/ref_CAG.fasta'

############################################################
#Basecall paths

bc_working_dir = '/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/XPCR/250504_B13_Basecall/BA-Basecall'
bc_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/raw/xPCR/240816_B13_pH_Rerun_2/20240816_1605_MN41475_AUB839_d3996fad/pod5'
bc_xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/raw/xPCR/240816_B13_pH_Rerun_2/reference/B13.fasta'
barcode_pair_csv = '/home/marchandlab/DataAnalysis/Kaplan/raw/xPCR/240816_B13_pH_Rerun_2/reference/DEMUX_B13.csv'
bc_model_file = '/home/marchandlab/github/kaplanna/xemora/models/240930_NTC_Models/GBC-BA-model_best.pt'

#bc_model_file = working_dir+'/model/model_best.pt'


############################################################
############################################################
train_model = False
basecall_reads = True
output_alignment_results = True
cutadapt_demux = True
raw_basecall_analysis = False
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
    cmd = f'python {demux_path} {bc_working_dir} {barcode_pair_csv}'
    os.system(cmd)    
    
#raw basecall analysis
if raw_basecall_analysis==True: 
    raw_bc_path = './lib/xr_raw_basecall_analysis.py'
    cmd = f'python {raw_bc_path} {bc_working_dir}'
    os.system(cmd)    
