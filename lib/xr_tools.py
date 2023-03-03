########################################################################
########################################################################
"""
xr_tools.py 

Title: Unpublished work

By: H. Kawabe, N. Kaplan, J. A. Marchand

Updated: 3/2/23
"""
########################################################################
########################################################################

import pandas as pd
import os
import glob
from pathlib import Path

def fetch_xna_pos(xm_header):
    pos=xm_header[xm_header.find('XPOS[')+5:-1].split('-')
    xpos = [x.split(':') for x in pos]
    return xpos


def xna_base_rc(xna_base, xna_bp): 
	for x in xna_bp: 
		if xna_base in x:
			if len(x)>1:
			    xx=list(x)
			    xx.remove(xna_base)
			    xrc=xx[0]
			else:
			   xrc=False
	return xrc


def check_xfasta_format(xfasta_file,standard_bases): 
    xfasta_header=False
    xna_in_sequence = False 
    with open(xfasta_file, "r") as fh:
        for line in fh:
            #Get header
            if line[0]=='>':
                header = line
                if 'XPOS[' in header: 
                    xfasta_header=True
                
            #Get sequence
            if line[0]!='>':
            
                #Make all upper case
                uline = line.upper()
                
                #Look for non-standard base
                diff = list(set(uline.replace('\n',''))-(set(standard_bases)))
                
                #This sequence contains an XNA
                if len(diff)>0:
                    xna_in_sequence = True 
    if xfasta_header==True and xna_in_sequence == False: 
        return True 
    else: 

        return False 



#Scans directory and subdirectory to get proper fast5 file path. Not explicitly handled with pod5 commands
def get_fast5_subdir(fast5_dir): 
    path = os.path.normpath(os.path.expanduser(fast5_dir))
    if os.path.exists(path):
        fast5_files = list(Path(path).rglob("*.fast5" ))
        if len(fast5_files)>0:
            fast5_subdir = os.path.dirname(fast5_files[0])
            print('Xemora [STATUS] - Found '+str(len(fast5_files))+' fast5 files in '+fast5_dir)
            return fast5_subdir
        else: 
            print('Xemora [ERROR] - Could not find Fast5 files in specified directory. Check .fast5 exist.')
            return False
    else: 
        print('Xemora [ERROR] - Could not find Fast5 directory. Check path')
        return False


#Check if working directory exists, if not create it. 
def check_make_dir(directory):
    directory = os.path.expanduser(directory)
    if not os.path.isdir(directory):
        os.makedirs(directory)
        print('Xemora [STATUS] - Required directory not found. Creating directory: '+directory)
    return directory


#Fast5 to pod5 conversion
def cod5_to_fast5(fast5_input, pod5_output):
    cmd = 'pod5 convert fast5 '+fast5_input+'/*.fast5 '+pod5_output
    os.system(cmd)
