import os
import glob
import subprocess
from collections import defaultdict

def filter_common_reads(demux_dir, output_common_reads):
    read_files = glob.glob(os.path.join(demux_dir, '*.fastq'))
    read_id_to_files = defaultdict(list)
    
    for read_file in read_files:
        with open(read_file, 'r') as f:
            for line in f:
                if line.startswith('@'):
                    read_id = line.split()[0]
                    read_id_to_files[read_id].append(read_file)
    
    with open(output_common_reads, 'w') as out:
        for read_id, files in read_id_to_files.items():
            if len(files) > 1:
                out.write(f"{read_id}\t{','.join(files)}\n")

def extract_read_ids(common_reads_file, output_read_ids):
    with open(common_reads_file, 'r') as infile, open(output_read_ids, 'w') as outfile:
        for line in infile:
            read_id = line.split()[0]
            outfile.write(f"{read_id}\n")

def remove_reads(demux_dir, read_ids_file):
    read_ids = set()
    with open(read_ids_file, 'r') as f:
        for line in f:
            read_ids.add(line.strip())
    
    read_files = glob.glob(os.path.join(demux_dir, '*.fastq'))
    for read_file in read_files:
        temp_file = read_file + '.tmp'
        with open(read_file, 'r') as infile, open(temp_file, 'w') as outfile:
            write_read = True
            for line in infile:
                if line.startswith('@'):
                    read_id = line.split()[0]
                    write_read = read_id not in read_ids
                if write_read:
                    outfile.write(line)
        os.replace(temp_file, read_file)

# Directories and file paths
demux_dir = "/home/xenolab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/240930_NTC_Phusion_Training_Testing/240930_NTC_Phusion_Withheld_Test_Set/demux"
output_common_reads = os.path.join(demux_dir, 'common_reads.txt')
output_read_ids = os.path.join(demux_dir, 'read_ids.txt')

# Step 1: Filter common reads
filter_common_reads(demux_dir, output_common_reads)
print(f"Common reads written to: {output_common_reads}")

# Step 2: Extract read IDs
extract_read_ids(output_common_reads, output_read_ids)
print(f"Read IDs written to: {output_read_ids}")

# Step 3: Remove reads from FASTQ files
remove_reads(demux_dir, output_read_ids)
print("Duplicate reads removed from FASTQ files.")

