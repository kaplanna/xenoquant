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
barcode_list = "./demux/barcode_files/NB_BARCODES.csv"


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

def load_full_barcode_list(barcode_list):
    """
    Load the full barcode list from a CSV file into a dictionary.
    
    Args:
        barcode_list (str): Path to the CSV file with BARCODE_NAME and BARCODE_SEQUENCE columns.
        
    Returns:
        dict: A dictionary where keys are BARCODE_NAME and values are BARCODE_SEQUENCE.
    """
    barcode_dict = {}
    with open(barcode_list, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            name = row['BARCODE_NAME']
            sequence = row['BARCODE_SEQUENCE']
            barcode_dict[name] = sequence
    return barcode_dict

def map_barcode_pairs(pair_csv, barcode_dict):
    """
    Map barcode pairs from names to sequences using the full barcode list dictionary.
    
    Args:
        pair_csv (str): Path to the CSV file with FWD_BARCODE and REV_BARCODE columns.
        barcode_dict (dict): Dictionary mapping BARCODE_NAME to BARCODE_SEQUENCE.
        
    Returns:
        dict: A dictionary where keys are pair names and values are tuples of forward and reverse sequences.
    """
    barcode_pairs = {}
    with open(pair_csv, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            fwd_name = row['FWD_BARCODE']
            rev_name = row['REV_BARCODE']
            pair_name = f"{fwd_name}_FWD_{rev_name}_REV"
            fwd_sequence = barcode_dict.get(fwd_name)
            rev_sequence = barcode_dict.get(rev_name)
            if not fwd_sequence or not rev_sequence:
                raise ValueError(f"Barcode name {fwd_name} or {rev_name} not found in the full barcode list.")
            barcode_pairs[pair_name] = (fwd_sequence, rev_sequence)
    return barcode_pairs


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

# Load the per-read modifications file
def load_per_read_modifications(file_path):
    return pd.read_csv(file_path, sep='\t')

# Merge the modifications data with the all_read_ids data
def merge_read_ids_with_modifications(modifications_df, read_ids_csv):
    # Load the read IDs CSV
    read_ids_df = pd.read_csv(read_ids_csv)
    
    # Merge on 'read_id', with left join to keep all the reads from the per-read modifications
    merged_df = pd.merge(modifications_df, read_ids_df, on='read_id', how='left')

    # If a read_id wasn't demultiplexed, barcode_pair will be NaN. Let's replace those with 'None'
    merged_df['barcode_pair'].fillna('None', inplace=True)

    return merged_df

# ave the merged dataframe to a new file
def save_merged_output(merged_df, output_file):
    merged_df.to_csv(output_file, sep='\t', index=False)

# Load the full barcode list
barcode_dict = load_full_barcode_list(barcode_list)

# Map barcode pairs to sequences
barcode_pairs = map_barcode_pairs(barcode_pair_csv, barcode_dict)

# Prepare to collect all read IDs
all_read_ids = []

# Process the BAM file with each barcode pair
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

# Create a function to aggregate results by barcode pair
def calculate_overall_results(merged_df):
    # Group by barcode_pair and count total alignments
    results = merged_df.groupby('barcode_pair').agg(
        Total_Reads=('read_id', 'size'),
        Number_of_1s=('class_pred', lambda x: (x == 1).sum()),
        Number_of_0s=('class_pred', lambda x: (x == 0).sum())
    ).reset_index()

    # Calculate the percentage of 1s and 0s
    results['Percentage_1'] = (results['Number_of_1s'] / results['Total_Reads']) * 100
    results['Percentage_0'] = (results['Number_of_0s'] / results['Total_Reads']) * 100

    return results

# Save the aggregated results to a CSV file
def save_overall_results(results_df, output_file):
    results_df.to_csv(output_file, index=False)


# Remove read IDs that appear in multiple barcode pairs and print the count
all_read_ids = remove_duplicate_read_ids(all_read_ids)

# Write all read IDs to a single CSV file
read_ids_csv = os.path.join(demux_dir, 'all_read_ids.csv')
with open(read_ids_csv, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["barcode_pair", "read_id"])
    writer.writerows(all_read_ids)


print('Xemora [Status]: Cutadapt processing, deduplication, read ID extraction, duplicate removal, and cleanup completed for all barcode pairs.')

# Load the per-read modifications file (assumes working_dir is already defined)
modifications_file = os.path.join(working_dir, 'remora_outputs', 'per-read_modifications.tsv')
modifications_df = load_per_read_modifications(modifications_file)

# Load the all_read_ids.csv generated by the previous script and merge it with the modifications
read_ids_csv = os.path.join(demux_dir, 'all_read_ids.csv')
merged_df = merge_read_ids_with_modifications(modifications_df, read_ids_csv)

# Save the result back to a new file
output_file = os.path.join(demux_dir, 'demux_per-read_modifications.tsv')
save_merged_output(merged_df, output_file)

print(f'Merged modifications and read IDs saved to {output_file}')

overall_results = calculate_overall_results(merged_df)

# Save the overall results to a CSV file
overall_results_file = os.path.join(demux_dir, 'overall_demux_results.csv')
save_overall_results(overall_results, overall_results_file)

print(f'Overall results saved to {overall_results_file}')
