import pysam
import pandas as pd
import logomaker
import matplotlib.pyplot as plt
from collections import Counter
import os
import sys
import subprocess
from tqdm import tqdm
from xr_tools  import *
from xr_params import *

print('Xemora [Status] - Initializing Xemora Raw Basecall Analysis.')

# Define working directories
working_dir = os.path.expanduser(sys.argv[1])
raw_basecall_output_dir = check_make_dir(os.path.join(working_dir, 'raw_basecall_analysis'))
raw_bc_output_file = os.path.join(raw_basecall_output_dir, 'raw_bc_output.csv')

demux_dir = check_make_dir(os.path.join(working_dir, 'demux'))
demux_read_IDs = os.path.join(demux_dir, 'all_read_ids.csv')

mod_dir = os.path.join(working_dir, 'preprocess')
mod_bam_dir = os.path.join(mod_dir, 'bam')
aligned_bam = os.path.join(mod_bam_dir, 'aligned.BAM')

ref_dir = os.path.join(working_dir, 'references')
bed_file = os.path.join(ref_dir, f'{mod_base}.bed')


def extract_positions_from_bed(bed_file):
    positions = []
    with open(bed_file, 'r') as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) >= 3:
                ref_name = cols[0]
                position = int(cols[1]) + 1  # BED is 0-based
                positions.append((ref_name, position))
    return positions


def reverse_complement(base):
    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "DEL": "DEL", "INS": "INS"}
    return complement.get(base, base)


def analyze_position(bam_file, reference_name, position):
    bam = pysam.AlignmentFile(bam_file, "rb")
    results = []
    target_ref_pos = position - 1

    for read in bam:
        if read.reference_name != reference_name or read.is_unmapped:
            continue

        read_id = read.query_name
        strand = '+' if not (read.flag & 16) else '-'

        ref_index = read.reference_start
        query_index = 0

        for cigar_op, length in read.cigartuples:
            if cigar_op == 0:
                for i in range(length):
                    current_ref_pos = ref_index + i
                    if current_ref_pos == target_ref_pos:
                        base = read.query_sequence[query_index + i]
                        if strand == '-':
                            base = reverse_complement(base)
                        results.append({'read_id': read_id, 'base': base, 'strand': strand})
                query_index += length
                ref_index += length
            elif cigar_op == 1:
                if ref_index == target_ref_pos:
                    results.append({'read_id': read_id, 'base': 'INS', 'strand': strand})
                query_index += length
            elif cigar_op == 2:
                if ref_index == target_ref_pos:
                    results.append({'read_id': read_id, 'base': 'DEL', 'strand': strand})
                ref_index += length
            elif cigar_op in [4, 5]:
                query_index += length

    bam.close()
    return results


def generate_logoplots_from_csv(csv_file):
    df = pd.read_csv(csv_file)
    all_base_cols = ['A', 'T', 'G', 'C', 'DEL', 'INS']
    nucleotide_cols = [col for col in all_base_cols if col in df.columns]

    df = df.reset_index()
    output_dir = os.path.dirname(csv_file)

    color_scheme = {
        'A': '#008000',  # green
        'T': '#ff0000',  # red
        'G': '#ffa500',  # orange
        'C': '#0000ff',  # blue
        'D': '#000000',  # black for DEL
        'I': '#808080'   # gray for INS
    }

    for strand in df['strand'].unique():
        df_strand = df[df['strand'] == strand].copy()
        freq_df = df_strand[nucleotide_cols].div(df_strand[nucleotide_cols].sum(axis=1), axis=0)
        freq_df = freq_df.reset_index(drop=True)
        freq_df.index = range(len(freq_df))

        freq_df.columns = [col if len(col) == 1 else col[0] for col in freq_df.columns]

        plt.figure(figsize=(12, 4))
        logo = logomaker.Logo(freq_df, shade_below=.5, fade_below=.5, font_name='Arial', color_scheme=color_scheme)
        logo.style_spines(visible=False)
        logo.style_spines(spines=['left', 'bottom'], visible=True)
        logo.ax.set_ylabel('Frequency')
        logo.ax.set_xlabel(f'Sample (Strand: {strand})')
        logo.ax.set_xticks(range(len(df_strand)))
        logo.ax.set_xticklabels(df_strand['sample_id'], rotation=45)

        plt.title(f'Sequence Logo - Strand: {strand}')
        plt.tight_layout()

        output_file = os.path.join(output_dir, f'seq_logo_strand_{strand}.pdf')
        plt.savefig(output_file)
        print(f"Saved logo plot to {output_file}")
        plt.close()


def main():
    positions = extract_positions_from_bed(bed_file)
    results = []
    for ref_name, position in tqdm(positions, desc="Processing Alignments", unit="alignments"):
        results.extend(analyze_position(aligned_bam, ref_name, position))

    df = pd.DataFrame(results)
    barcode_df = pd.read_csv(demux_read_IDs, sep=',')

    required_columns = {"sample_id", "barcode_pair", "read_id"}
    if not required_columns.issubset(set(barcode_df.columns)):
        print(f"Error: all_read_ids.csv must contain columns {required_columns}, but found: {barcode_df.columns}")
        return

    df = df.merge(barcode_df, on='read_id', how='left')

    summary_df = df.groupby(['sample_id', 'barcode_pair', 'strand'])['base'].value_counts().unstack(fill_value=0)

    desired_order = ['INS', 'DEL', 'A', 'T', 'G', 'C']
    summary_df = summary_df.reindex(columns=[col for col in desired_order if col in summary_df.columns], fill_value=0)

    summary_df.to_csv(raw_bc_output_file)
    print(f"Summary results saved to {raw_bc_output_file}")

    generate_logoplots_from_csv(raw_bc_output_file)


if __name__ == "__main__":
    main()

