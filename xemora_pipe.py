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

working_dir = '/home/marchandlab/DataAnalysis/Kaplan/training/10.4.1/230928_GAP_Train/ZC_Train_extratest'
xna_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/fast5/10.4.1/230928_GAP_Train/GAP_basecall_50fast5'
xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/training/10.4.1/reference/GAP_allseq_wadapt.fa'
dna_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/fast5/10.4.1/230928_GAP_Train/GAP_basecall_50fast5'
dna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/training/10.4.1/reference/GAP_allseq_wadapt.fa'
############################################################

############################################################
#Basecall paths

bc_working_dir = '/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/230928_GAP_Train_Test/ZC_Basecall'
bc_fast5_dir = '/home/marchandlab/DataAnalysis/Kaplan/fast5/10.4.1/230928_GAP_Train/GAP_basecall_92fast5'
bc_xna_ref_fasta = '/home/marchandlab/DataAnalysis/Kaplan/training/10.4.1/reference/GAP_allseq_wadapt.fa'
bc_model_file = '/home/marchandlab/DataAnalysis/Kaplan/training/10.4.1/230928_GAP_Train/ZC_Train/model/model_best.pt'

############################################################

############################################################
train_model =  True
basecall_reads =  False
output_basecall_results = False
############################################################

#conda activate nanopore-re3


        
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
    
