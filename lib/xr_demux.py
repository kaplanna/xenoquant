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
barcode_pair_csv = os.path.expanduser(sys.argv[2])
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
    result = subprocess.run(cmd, check=True, text=True, capture_output=True)
    print(result.stderr)  # short summary


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
    Also, extract the sample ID.
    
    Args:
        pair_csv (str): Path to the CSV file with SAMPLE_ID, FWD_BARCODE, and REV_BARCODE columns.
        barcode_dict (dict): Dictionary mapping BARCODE_NAME to BARCODE_SEQUENCE.
        
    Returns:
        dict: A dictionary where keys are pair names, values are tuples of (sample ID, forward sequence, reverse sequence).
    """
    barcode_pairs = {}
    with open(pair_csv, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            sample_id = row['SAMPLE_ID']
            fwd_name = row['FWD_BARCODE']
            rev_name = row['REV_BARCODE']
            pair_name = f"{fwd_name}_FWD_{rev_name}_REV"
            fwd_sequence = barcode_dict.get(fwd_name)
            rev_sequence = barcode_dict.get(rev_name)
            if not fwd_sequence or not rev_sequence:
                raise ValueError(f"Barcode name {fwd_name} or {rev_name} not found in the full barcode list.")
            
            barcode_pairs[pair_name] = (sample_id, fwd_sequence, rev_sequence)
    return barcode_pairs



def extract_read_ids(fastq_file, barcode_pair, sample_id, read_ids_list):
    with open(fastq_file, 'r') as infile:
        for line in infile:
            if line.startswith('@'):
                read_id = line.split()[0][1:]  # Extract the read ID
                read_ids_list.append([sample_id, barcode_pair, read_id])


def remove_duplicate_read_ids(all_read_ids):
    read_id_counts = Counter(read_id for _, _, read_id in all_read_ids)
    duplicates = [entry for entry in all_read_ids if read_id_counts[entry[2]] > 1]
    unique_entries = [entry for entry in all_read_ids if read_id_counts[entry[2]] == 1]
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
    
def calculate_overall_results(merged_df, barcode_pairs):
    """
    Calculate overall results grouped by sample ID and barcode pair.
    
    Args:
        merged_df (DataFrame): Merged DataFrame with barcode pair and read classifications.
        barcode_pairs (dict): Dictionary with barcode pairs and their corresponding sample IDs.
        
    Returns:
        DataFrame: Aggregated results including sample ID, barcode pair, and classification counts.
    """
    # Add sample ID column by mapping from barcode_pairs dictionary
    merged_df['sample_id'] = merged_df['barcode_pair'].map(
        lambda x: barcode_pairs.get(x, ('Unknown',))[0]
    )

    # Group by sample ID and barcode pair
    results = merged_df.groupby(['sample_id', 'barcode_pair']).agg(
        Total_Reads=('read_id', 'size'),
        Number_of_1s=('class_pred', lambda x: (x == 1).sum()),
        Number_of_0s=('class_pred', lambda x: (x == 0).sum())
    ).reset_index()

    # Calculate percentages
    results['Percentage_1'] = (results['Number_of_1s'] / results['Total_Reads']) * 100
    results['Percentage_0'] = (results['Number_of_0s'] / results['Total_Reads']) * 100

    return results

def save_overall_results(results_df, output_file):
    results_df.to_csv(output_file, index=False)
    
def extract_p_class1(prob_string):
    """Extract P(class1) from 'p0,p1'"""
    try:
        return float(prob_string.split(",")[1])
    except Exception:
        return np.nan


def apply_hard_decision_threshold(df, threshold):
    """
    Reassign basecalls using a hard decision threshold (>0.5).
    No reads are dropped.
    """
    df = df.copy()
    df["p_class1"] = df["class_probs"].apply(extract_p_class1)
    df["class_pred_reassigned"] = (df["p_class1"] >= threshold).astype(int)
    return df


def calculate_overall_results_reassigned(df):
    """
    Aggregate results using reassigned basecalls
    """
    results = (
        df.groupby(["sample_id", "barcode_pair"])
          .agg(
              Total_Reads=("read_id", "size"),
              Number_of_1s=("class_pred_reassigned", lambda x: (x == 1).sum()),
              Number_of_0s=("class_pred_reassigned", lambda x: (x == 0).sum())
          )
          .reset_index()
    )

    results["Percentage_1"] = 100 * results["Number_of_1s"] / results["Total_Reads"]
    results["Percentage_0"] = 100 * results["Number_of_0s"] / results["Total_Reads"]

    return results

    
    
# Load the full barcode list
barcode_dict = load_full_barcode_list(barcode_list)

# Map barcode pairs to sequences
barcode_pairs = map_barcode_pairs(barcode_pair_csv, barcode_dict)

# Prepare to collect all read IDs
all_read_ids = []

if RERUN_DEMUX:
    # Process the BAM file with each barcode pair
    for barcode_pair, (sample_id, barcode1, barcode2) in barcode_pairs.items():
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

        dedup_output = os.path.join(fastq_dir, f'{barcode_pair}.fastq')
        print(f'Combining and deduplicating output for barcode pair: {barcode_pair}')
        combine_and_deduplicate_fastq(output_file1, output_file2, dedup_output)

        extract_read_ids(dedup_output, barcode_pair, sample_id, all_read_ids)

        os.remove(output_file1)
        os.remove(output_file2)

    # Remove duplicates and save all_read_ids.csv
    all_read_ids = remove_duplicate_read_ids(all_read_ids)
    read_ids_csv = os.path.join(demux_dir, 'all_read_ids.csv')

    if not all_read_ids:
        print("Warning: No read IDs collected, skipping write to all_read_ids.csv")
    else:
        with open(read_ids_csv, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["sample_id", "barcode_pair", "read_id"])
            writer.writerows(all_read_ids)

    print('Xemora [Status]: Cutadapt demux completed.')

else:
    # Skip demux, just load existing read_ids
    read_ids_csv = os.path.join(demux_dir, 'all_read_ids.csv')
    if not os.path.exists(read_ids_csv):
        raise FileNotFoundError(
            f"No existing all_read_ids.csv found at {read_ids_csv}. "
            "Set RERUN_DEMUX=True to generate it."
        )
    print(f'Xemora [Status]: Skipping demux. Using existing {read_ids_csv}.')

# -------------------------------------------------------------------
# Analysis section (unchanged)
# -------------------------------------------------------------------

# Load the per-read modifications file (assumes working_dir is already defined)
modifications_file = os.path.join(working_dir, 'remora_outputs', 'per-read_modifications.tsv')
modifications_df = load_per_read_modifications(modifications_file)

# Merge with read IDs
merged_df = merge_read_ids_with_modifications(modifications_df, read_ids_csv)

output_file = os.path.join(demux_dir, 'demux_per-read_modifications.tsv')

# Keep only the desired columns
keep_cols = [
    "read_id",
    "read_focus_base",
    "label",
    "class_pred",
    "class_probs",
    "ref_start_pos",
    "ref_length",
    "basecalled_sequence",
    "sample_id",
    "barcode_pair"
]

slim_df = merged_df[keep_cols].copy()
slim_df.to_csv(output_file, sep='\t', index=False)

print(f'Merged modifications (slimmed) saved to {output_file}')


# Aggregate overall results
overall_results = calculate_overall_results(merged_df, barcode_pairs)

# Save overall results
overall_results_file = os.path.join(demux_dir, 'overall_demux_results.csv')
save_overall_results(overall_results, overall_results_file)
print(f'Overall results saved to {overall_results_file}')

# -------------------------------------------------------------------
# Optional: Decision threshold reassignment
# -------------------------------------------------------------------

if USE_DECISION_THRESHOLD:
    print(
        f"Xemora [Status]: Reassigning basecalls with decision threshold = {DECISION_THRESHOLD}"
    )

    reassigned_df = apply_hard_decision_threshold(
        slim_df,
        DECISION_THRESHOLD
    )

    # ------------------------------------------------------------
    # Save per-read reassigned output (NEW FILE)
    # ------------------------------------------------------------
    reassigned_per_read = os.path.join(
        demux_dir,
        f"demux_per-read_modifications_recalled_T{DECISION_THRESHOLD:.2f}.tsv"
    )

    reassigned_df.to_csv(reassigned_per_read, sep="\t", index=False)

    print(f"Reassigned per-read file saved to {reassigned_per_read}")

    # ------------------------------------------------------------
    # Recompute overall results (NEW FILE)
    # ------------------------------------------------------------
    reassigned_overall = calculate_overall_results_reassigned(reassigned_df)

    reassigned_overall_file = os.path.join(
        demux_dir,
        f"overall_demux_results_recalled_T{DECISION_THRESHOLD:.2f}.csv"
    )

    reassigned_overall.to_csv(reassigned_overall_file, index=False)

    print(f"Reassigned overall results saved to {reassigned_overall_file}")

    # ------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------
    changed = (
        reassigned_df["class_pred"] !=
        reassigned_df["class_pred_reassigned"]
    ).sum()

    print(
        f"Decision threshold changed {changed} / {len(reassigned_df)} reads "
        f"relative to Remora argmax."
    )


