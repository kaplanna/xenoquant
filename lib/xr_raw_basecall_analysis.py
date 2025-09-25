import pysam
import pandas as pd
import logomaker
import matplotlib.pyplot as plt
import numpy as np
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
modification_filter_file = os.path.join(working_dir, 'remora_outputs', 'per-read_modifications.tsv') \
    if FILTER_BY_CLASS_0 else None



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


def generate_logoplots_from_csv(csv_file, suffix=""):
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

        output_file = os.path.join(output_dir, f'seq_logo_strand_{strand}{suffix}.pdf')
        plt.savefig(output_file)
        print(f"Saved logo plot to {output_file}")
        plt.close()



def plot_tornado(df, output_file):
    """
    Create a tornado plot where Call=0 fractions go left and Call=1 fractions go right.
    Plot order matches the order of rows in df (overall_summary).
    """
    desired_order = ['A', 'T', 'G', 'C', 'INS', 'DEL']
    color_scheme = {
        'A': '#008000',  # green
        'T': '#ff0000',  # red
        'G': '#ffa500',  # orange
        'C': '#0000ff',  # blue
        'INS': '#808080',# gray
        'DEL': '#000000' # black
    }

    # Extract unique (sample_id, strand) pairs in the same order they appear in df
    order_keys = df[['sample_id','strand']].drop_duplicates().values.tolist()

    fig, ax = plt.subplots(figsize=(10, 6))
    y_labels = []
    y_pos = np.arange(len(order_keys))

    for i, (sample_id, strand) in enumerate(order_keys):
        y_labels.append(f"{sample_id}|{strand}")
        group = df[(df['sample_id']==sample_id) & (df['strand']==strand)]

        row0 = group[group['Call'] == 0]
        row1 = group[group['Call'] == 1]

        left_start = 0
        right_start = 0

        for base in desired_order:
            col = f"weighted_{base}"
            if col not in df.columns:
                continue

            val0 = row0[col].values[0] if not row0.empty else 0
            val1 = row1[col].values[0] if not row1.empty else 0

            # Call 0 → left (negative)
            ax.barh(i, -val0, left=left_start, color=color_scheme[base], label=base)
            left_start -= val0

            # Call 1 → right (positive)
            ax.barh(i, val1, left=right_start, color=color_scheme[base], label=base)
            right_start += val1

    # Set y-axis labels in the same order as df
    ax.set_yticks(y_pos)
    ax.set_yticklabels(y_labels, fontsize=8)
    ax.set_xlabel("Weighted relative frequency")
    ax.set_title("Tornado Plot of Base Frequencies (Call 0 vs 1)")
    ax.axvline(0, color='black', linewidth=0.8)

    # Deduplicate legend
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), title="Base", bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    print(f"[Status] Tornado plot saved to {output_file}")




def main():
    positions = extract_positions_from_bed(bed_file)
    combined_summaries = []

    for class_val in [0, 1]:
        print(f"\n[Status] Processing for class_pred == {class_val}")
        results = []
        for ref_name, position in tqdm(
            positions,
            desc=f"Alignments for class_pred == {class_val}",
            unit="alignments"
        ):
            results.extend(analyze_position(aligned_bam, ref_name, position))

        df = pd.DataFrame(results)

        # Filter based on class_pred
        if os.path.isfile(modification_filter_file):
            try:
                mod_df = pd.read_csv(
                    modification_filter_file,
                    sep='\t',
                    usecols=['read_id', 'class_pred']
                )
                filtered_ids = set(mod_df[mod_df['class_pred'] == class_val]['read_id'])
                df = df[df['read_id'].isin(filtered_ids)]
            except Exception as e:
                print(f"Error reading or processing modification filter file: {e}")
                continue

        # Merge with demux info
        barcode_df = pd.read_csv(demux_read_IDs, sep=',')
        required_columns = {"sample_id", "barcode_pair", "read_id"}
        if not required_columns.issubset(set(barcode_df.columns)):
            print(f"Error: all_read_ids.csv must contain columns {required_columns}, but found: {barcode_df.columns}")
            continue

        df = df.merge(barcode_df, on='read_id', how='left')

        # ---- summary table per alignment/demux ----
        summary_df = (
            df.groupby(['sample_id', 'barcode_pair', 'strand'])['base']
              .value_counts()
              .unstack(fill_value=0)
        )

        desired_order = ['INS', 'DEL', 'A', 'T', 'G', 'C']
        summary_df = summary_df.reindex(
            columns=[col for col in desired_order if col in summary_df.columns],
            fill_value=0
        )

        # Add Call column
        summary_df['Call'] = class_val

        # Add Total + relative frequencies
        summary_df['Total'] = summary_df[desired_order].sum(axis=1)
        for col in desired_order:
            if col in summary_df.columns:
                summary_df[f'rel_{col}'] = summary_df[col] / summary_df['Total']

        # Save per-class summary
        class_output_file = os.path.join(
            raw_basecall_output_dir,
            f'raw_bc_output_{class_val}.csv'
        )
        summary_df.to_csv(class_output_file)
        print(f"Summary results saved to {class_output_file}")

        generate_logoplots_from_csv(class_output_file, suffix=f"_{class_val}")

        # Keep for overall summary
        combined_summaries.append(summary_df)

        # ---- AFTER loop: write combined summary + tornado plot ----
    if combined_summaries:
        overall_summary = pd.concat(combined_summaries).reset_index()

        # Compute weights per (sample_id, barcode_pair, strand)
        desired_order = ['INS', 'DEL', 'A', 'T', 'G', 'C']
        def apply_weighting(group):
            total_all = group['Total'].sum()
            if total_all == 0:
                return group
            weight = group['Total'] / total_all
            for base in desired_order:
                rel_col = f'rel_{base}'
                if rel_col in group.columns:
                    group[f'weighted_{base}'] = group[rel_col] * weight
            group['Call_weight'] = weight
            return group

        overall_summary = overall_summary.groupby(
            ['sample_id','barcode_pair','strand'], group_keys=False
        ).apply(apply_weighting)

        overall_output_file = os.path.join(
            raw_basecall_output_dir,
            "raw_bc_output_overall.csv"
        )
        overall_summary.to_csv(overall_output_file, index=False)
        print(f"[Status] Overall summary with weighted frequencies saved to {overall_output_file}")

        tornado_output = os.path.join(
            raw_basecall_output_dir,
            "tornado_plot.pdf"
        )
        plot_tornado(overall_summary, tornado_output)






if __name__ == "__main__":
    main()

