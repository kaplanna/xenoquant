import pandas as pd
import numpy as np
import os
import sys
import subprocess
import csv
from collections import Counter
from xr_tools  import *
from xr_params import *

print('Xemora [Status] - Initializing Xemora Demux.')
working_dir = os.path.expanduser(sys.argv[1])
demux_dir = check_make_dir(os.path.join(working_dir, 'demux'))
fastq_dir = check_make_dir(os.path.join(demux_dir, 'fastq'))
mod_dir = os.path.join(working_dir,'preprocess')
mod_bam_dir = os.path.join(mod_dir,'bam')
bc_bam = os.path.join(mod_bam_dir,'bc.bam')
barcode_file = "./demux/xpcr_barcodes.fasta"

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in reversed(seq)])

def run_cutadapt(input_fastq, adapter, output_file, error_rate, min_overlap, min_len, max_len):
    cmd = [
        'cutadapt',
        '-g', adapter,  # Forward barcode
        '-m', str(min_len),  # Minimum length
        '-M', str(max_len),  # Maximum length
        '--no-trim',  # Do not trim off the barcode sequence
        '--discard-untrimmed',  # Discard untrimmed reads
        '-e', str(error_rate),  # Error rate
        '-O', str(min_overlap),  # Minimum overlap
        '-o', output_file,  # Output file
        input_fastq  # Input FASTQ file
    ]
    subprocess.run(cmd, check=True)

def combine_and_deduplicate_fastq(file1, file2, output_file):
    unique_reads = set()
    with open(file1, 'r') as infile1, open(file2, 'r') as infile2, open(output_file, 'w') as outfile:
        for infile in [infile1, infile2]:
            for title, seq, plus, qual in zip(*[infile]*4):
                if seq not in unique_reads:
                    unique_reads.add(seq)
                    outfile.writelines([title, seq, plus, qual])

def parse_barcodes(barcode_file):
    forward_barcodes = {}
    reverse_barcodes = {}
    with open(barcode_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                name = line.strip()[1:]
                sequence = next(f).strip()
                if 'FWD' in name:
                    forward_barcodes[name] = sequence
                elif 'REV' in name:
                    reverse_barcodes[name] = sequence
    return forward_barcodes, reverse_barcodes

def generate_barcode_pairs(forward_barcodes, reverse_barcodes):
    pairs = {}
    for fwd_name, fwd_seq in forward_barcodes.items():
        for rev_name, rev_seq in reverse_barcodes.items():
            pair_name = f"{fwd_name}_{rev_name}"
            pairs[pair_name] = (fwd_seq, rev_seq)
    return pairs



def extract_read_ids(fastq_file, barcode_pair, read_ids_list):
    with open(fastq_file, 'r') as infile:
        for line in infile:
            if line.startswith('@'):
                read_id = line.split()[0][1:]  # Extract the read ID
                read_ids_list.append([barcode_pair, read_id])

def remove_duplicate_read_ids(all_read_ids):
    read_id_counts = Counter(read_id for _, read_id in all_read_ids)
    duplicates = [entry for entry in all_read_ids if read_id_counts[entry[1]] > 1]
    unique_entries = [entry for entry in all_read_ids if read_id_counts[entry[1]] == 1]
    print(f"Number of read IDs found in multiple files: {len(duplicates)}")
    return unique_entries


error_rate = 0.20
min_overlap = 14
min_len = 120
max_len = 200


# Parse the barcodes from the file
forward_barcodes, reverse_barcodes = parse_barcodes(barcode_file)

# Generate all forward-reverse barcode pairs
barcode_pairs = generate_barcode_pairs(forward_barcodes, reverse_barcodes)

# Prepare to collect all read IDs
all_read_ids = []

# Process the bam file with each barcode pair
for barcode_pair, (barcode1, barcode2) in barcode_pairs.items():
    barcode2_rc = reverse_complement(barcode2)
    barcode1_rc = reverse_complement(barcode1)
    combined_adapter1 = f'{barcode1}...{barcode2_rc}'
    combined_adapter2 = f'{barcode2}...{barcode1_rc}'

    output_file1 = os.path.join(fastq_dir, f'{barcode_pair}_1.fastq')
    output_file2 = os.path.join(fastq_dir, f'{barcode_pair}_2.fastq')

    print(f'Running cutadapt for barcode pair: {barcode_pair}, adapter 1')
    run_cutadapt(bc_bam, combined_adapter1, output_file1, error_rate, min_overlap, min_len, max_len)

    print(f'Running cutadapt for barcode pair: {barcode_pair}, adapter 2')
    run_cutadapt(bc_bam, combined_adapter2, output_file2, error_rate, min_overlap, min_len, max_len)

    # Combine and deduplicate outputs
    dedup_output = os.path.join(fastq_dir, f'{barcode_pair}.fastq')
    print(f'Combining and deduplicating output for barcode pair: {barcode_pair}')
    combine_and_deduplicate_fastq(output_file1, output_file2, dedup_output)

    # Extract read IDs and append to the list
    extract_read_ids(dedup_output, barcode_pair, all_read_ids)

    # Delete the non-deduplicated files
    os.remove(output_file1)
    os.remove(output_file2)

# Remove read IDs that appear in multiple barcode pairs and print the count
all_read_ids = remove_duplicate_read_ids(all_read_ids)

# Write all read IDs to a single CSV file
read_ids_csv = os.path.join(demux_dir, 'all_read_ids.csv')
with open(read_ids_csv, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["barcode_pair", "read_id"])
    writer.writerows(all_read_ids)



print('Xemora [Status]: Cutadapt processing, deduplication, read ID extraction, duplicate removal, and cleanup completed for all barcode pairs.')
