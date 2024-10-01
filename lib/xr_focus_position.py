########################################################################
########################################################################
"""
xr_focus_position.py 

Title: Unpublished work

xr_focus_position.py takes outputs from realign and calculates per position 
averages for plotting

By: J. Sumabat, G. Loia, J. A. Marchand

Updated: 08/22/24
"""
########################################################################
########################################################################

########################################################################
# Imports
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sys
import re
import ast  # For safely evaluating strings as expressions
from xr_params import *
########################################################################

def load_data(input_file):
    """
    Load the TSV file into a DataFrame.

    Parameters:
    - input_file: str, the path to the input TSV file

    Returns:
    - df: pandas DataFrame, the loaded data
    """
    return pd.read_csv(input_file, sep='\t')

def filter_data(df, prob_filter):
    """
    Apply a probability filter to the DataFrame.

    Parameters:
    - df: pandas DataFrame, the input data
    - prob_filter: float, the threshold for filtering prob_DNA values

    Returns:
    - df_filtered: pandas DataFrame, the filtered data
    """
    return df[(df['prob_DNA'] >= prob_filter) | (df['prob_DNA'] <= 1 - prob_filter)]

def parse_cigar(cigar):
    """
    Parse the CIGAR string and return a list of operations.
    Each operation is a tuple of (operation, length).

    Parameters:
    - cigar: str, the CIGAR string

    Returns:
    - List of tuples: [(operation, length), ...]
    """
    return [(op[-1], int(op[:-1])) for op in re.findall(r'(\d+[MIDNSHP=X])', cigar)]

def convert_q_score_to_list(q_score_str):
    """
    Convert a q_score string in the format "array('B', [8, 10, 10, ...])"
    into a list of integers.

    Parameters:
    - q_score_str: str, the string representation of q_score

    Returns:
    - List[int]: A list of integers representing the q_scores
    """
    try:
        match = re.search(r"array\('B', \[(.*?)\]\)", q_score_str)
        if match:
            q_score_list = ast.literal_eval(f"[{match.group(1)}]")
            return q_score_list
        else:
            return None
    except Exception as e:
        print(f"Error converting q_score: {e}")
        return None

def resolve_base_position(cigar, ref_start_pos, read_focus_base, basecalled_sequence, q_scores, strand_flag, ref_length):
    """
    Resolve the base and corresponding quality score at a given focus position 
    by interpreting the CIGAR string.

    Parameters:
    - cigar: str, the CIGAR string
    - ref_start_pos: int, the start position on the reference (0-based)
    - read_focus_base: int, the target base on the reference (0-based)
    - basecalled_sequence: str, the basecalled sequence
    - q_scores: list or array of quality scores corresponding to the basecalled sequence
    - strand_flag: int, the SAM flag indicating the strand (0 = forward, 16 = reverse)

    Returns:
    - tuple: (char, int), where 'char' is the base called at the focus position,
             and 'int' is the corresponding quality score.
             Returns (None, None) if the base cannot be resolved.
    """
    
    def complement_base(base):
        """Return the complement of a base (for reverse strand reads)."""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return complement.get(base, base)

    operations = parse_cigar(cigar)
    sequence_index = 0
    reference_index = ref_start_pos

    if strand_flag == 16:
        target_position = ref_length - read_focus_base -1
    else:
        target_position = read_focus_base

    # Iterate through CIGAR operations to find the match with respect to the reference
    for op, length in operations:
        if op == 'S':  # Soft clip
            sequence_index += length  # Skip these bases in the sequence
        elif op == 'M':  # Match/mismatch
            if reference_index <= target_position < reference_index + length:
                # Calculate the offset within the sequence
                offset = target_position - reference_index
                base = basecalled_sequence[sequence_index + offset]
                q_score = q_scores[sequence_index + offset]  # Get the quality score
                
                # If the read is on the reverse strand, complement the base
               
                if strand_flag == 16:
                    base = complement_base(base)
               
                return base, q_score

            # Move along the sequence and reference
            reference_index += length
            sequence_index += length
        elif op == 'I':  # Insertion
            # Insertion doesn't affect reference position, just skip sequence
            sequence_index += length
        elif op == 'D' or op == 'N':  # Deletion/skip
            # Deletion affects reference but not the sequence
            reference_index += length
        elif op == 'H' or op == 'P':  # Hard clip or padding, just skip
            continue

    return None, None


def calculate_base_fractions(df):
    """
    Calculate the fraction of each base called at each reference position.

    Parameters:
    - df: pandas DataFrame, the filtered data with 'called_base' column

    Returns:
    - base_fractions: pandas DataFrame, where each row is a reference position,
      and columns are fractions of bases called (A, T, C, G, etc.)
    """
    df_filtered = df.dropna(subset=['called_base'])

    base_counts = df_filtered.groupby(['normalized_read_focus_base', 'called_base']).size().unstack(fill_value=0)
    total_counts = base_counts.sum(axis=1)
    base_fractions = base_counts.div(total_counts, axis=0)

    return base_fractions

def calculate_shannon_entropy(base_fractions):
    """
    Calculate the Shannon entropy for each reference position based on the base fractions.

    Parameters:
    - base_fractions: pandas DataFrame, columns are base fractions (A, T, C, G)

    Returns:
    - entropy: pandas Series, the Shannon entropy for each reference position
    """
    base_columns = ['A', 'C', 'G', 'T']
    base_fractions = base_fractions.reindex(columns=base_columns, fill_value=0)

    def entropy_formula(row):
        entropy = -sum(p * np.log2(p) for p in row if p > 0)
        return entropy

    entropy = base_fractions.apply(entropy_formula, axis=1)

    return entropy

def calculate_summary_statistics(df):
    """
    Calculate various summary statistics and combine them into a single DataFrame.

    This function performs several group-level aggregations on a filtered dataset, including:
    - Value counts for class predictions.
    - Mean and standard deviation for probability columns.
    - Odds ratio and log10 odds ratio.
    - Average and standard deviation of q-scores.
    - Shannon entropy for base fractions.

    Parameters:
    - df: pandas DataFrame
        The input DataFrame containing read information, including columns for 'prob_DNA', 'prob_XNA', 
        'q_score', 'class_pred', 'normalized_read_focus_base', and other relevant data.

    Returns:
    - summary: pandas DataFrame
        The resulting DataFrame that contains the summary statistics for each 'normalized_read_focus_base'.
    """
    
    # Group by 'normalized_read_focus_base' and calculate normalized value counts for 'class_pred'
    summary = df.groupby('normalized_read_focus_base')['class_pred'].value_counts(normalize=True).unstack(fill_value=0)

    # Calculate the mean and standard deviation for 'prob_DNA' and 'prob_XNA'
    avg_std = df.groupby('normalized_read_focus_base').agg({
        'prob_DNA': ['mean', 'std'],  # DNA probability
        'prob_XNA': ['mean', 'std']   # XNA probability
    })
    
    # Flatten the multi-level columns resulting from aggregation
    avg_std.columns = ['_'.join(col).strip() for col in avg_std.columns.values]

    # Calculate odds ratio and log10 odds ratio
    df['odds_ratio'] = df['prob_XNA'] / df['prob_DNA']  # Ratio of XNA to DNA probability
    df['log10_odds_ratio'] = np.log10(df['odds_ratio'])  # Log10 transformation for better interpretability

    # Group by 'normalized_read_focus_base' and calculate mean and std for odds ratio and log10 odds ratio
    odds_ratio_avg_std = df.groupby('normalized_read_focus_base')['odds_ratio'].agg(['mean', 'std'])
    log10_odds_ratio_avg_std = df.groupby('normalized_read_focus_base')['log10_odds_ratio'].agg(['mean', 'std'])

    # Calculate a "first odds ratio" using the mean values of prob_XNA and prob_DNA
    df['avg_first_odds_ratio'] = avg_std['prob_XNA_mean'] / avg_std['prob_DNA_mean']
    
    # Create a DataFrame to store the average and (placeholder) std of first odds ratio
    avg_odds_ratio_avg_std = pd.DataFrame(columns=['mean', 'std'])
    avg_odds_ratio_avg_std['mean'] = df['avg_first_odds_ratio']
    avg_odds_ratio_avg_std['std'] = np.nan  # No standard deviation for first odds ratio at this step

    # Convert the q-scores from a string format to list format for processing
    df['q_score_list'] = df['q_score'].apply(lambda x: convert_q_score_to_list(x) if isinstance(x, str) else x)

    # Resolve the base and q-score using the custom resolve_base_position function
    df['called_base'], df['q_score'] = zip(*df.apply(
        lambda row: resolve_base_position(
            row['cigar_string'], row['ref_start_pos'], row['read_focus_base'], 
            row['basecalled_sequence'], row['q_score_list'], row['flag'], row['ref_length']
        ), axis=1
    ))

    # Convert the q_score column to numeric, coercing any invalid values to NaN
    df['q_score'] = pd.to_numeric(df['q_score'], errors='coerce')

    # Calculate base fractions (A, C, G, T) using a custom function
    base_fractions = calculate_base_fractions(df)

    # Calculate Shannon entropy from the base fractions (measure of uncertainty in base distribution)
    entropy = calculate_shannon_entropy(base_fractions)

    # Calculate the average and standard deviation of q-scores for each group
    avg_q_score = df.groupby('normalized_read_focus_base')['q_score'].mean()
    std_q_score = df.groupby('normalized_read_focus_base')['q_score'].std()

    # Join calculated statistics into the summary DataFrame
    summary = summary.join(avg_std)  # Join avg and std for DNA/XNA probabilities
    summary = summary.join(odds_ratio_avg_std.rename(columns={'mean': 'avg_odds_ratio', 'std': 'std_odds_ratio'}))  # Odds ratio stats
    summary = summary.join(log10_odds_ratio_avg_std.rename(columns={'mean': 'avg_log10_odds_ratio', 'std': 'std_log10_odds_ratio'}))  # Log10 odds ratio stats
    summary = summary.join(avg_odds_ratio_avg_std.rename(columns={'mean': 'avg_first_odds_ratio', 'std': 'avg_first_std_odds_ratio'}))  # First odds ratio
    summary = summary.join(base_fractions)  # Base fractions (A, C, G, T)
    summary['shannon_entropy'] = entropy  # Shannon entropy
    summary['avg_q_score'] = avg_q_score  # Average q-score
    summary['std_q_score'] = std_q_score  # Standard deviation of q-score

    # Rename columns for clarity and consistency
    summary.columns = [
        'class_pred_0', 'class_pred_1', 'avg_prob_DNA', 'std_prob_DNA', 
        'avg_prob_XNA', 'std_prob_XNA', 'avg_odds_ratio', 'std_odds_ratio', 
        'avg_log10_odds_ratio', 'std_log10_odds_ratio', 'avg_first_odds_ratio', 
        'avg_first_std_odds_ratio', 'A', 'C', 'G', 'T', 'shannon_entropy', 
        'avg_q_score', 'std_q_score'  # Include standard deviation of q-scores
    ]

    return summary


def save_summary_to_csv(summary, output_file):
    """
    Save the summary DataFrame to a TSV file.

    Parameters:
    - summary: pandas DataFrame, the summary statistics
    - output_file: str, the path to save the TSV file

    Returns:
    - None
    """
    summary.to_csv(output_file, sep='\t')
    print(f"Xemora [STATUS] - Position analysis complete. Per sequence results saved to {output_file}")

def calculate_shannon_entropy_for_error_rate(correct_p=0.85, error_p=0.05):
    """
    Calculate the Shannon entropy based on the correct call probability and error call probabilities.
    
    Parameters:
    - correct_p: float, probability of the correct base being called (e.g., 0.85 for 15% error rate)
    - error_p: float, probability of an incorrect base being called (e.g., 0.05 for 15% error rate)
    
    Returns:
    - entropy: float, calculated Shannon entropy
    """
    entropy = -(correct_p * np.log2(correct_p) + 3 * error_p * np.log2(error_p))
    return entropy

def plot_summary(summary, output_graph_file_path, correct_p=0.85, error_p=0.05, display_vis=False):
    """
    Plot the summary statistics using matplotlib and add a red horizontal line
    for the Shannon entropy plot.

    Parameters:
    - summary: pandas DataFrame, the summary statistics
    - output_graph_file_path: str, path to save the output plot
    - correct_p: float, probability of the correct base being called (e.g., 0.85 for 15% error rate)
    - error_p: float, probability of an incorrect base being called (e.g., 0.05 for 15% error rate)
    - display_vis: bool, whether to display the plot in the notebook or script (default: False)
    
    Returns:
    - None
    """
    # Specify the columns to be plotted
    plot_columns = [
        'class_pred_0', 'class_pred_1', 'avg_prob_DNA', 'avg_prob_XNA', 
        'avg_odds_ratio', 'avg_log10_odds_ratio', 'avg_first_odds_ratio',
        'shannon_entropy', 'avg_q_score'
    ]

    # Set up the subplot grid (3 rows, 3 columns) with a larger figure size
    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(15, 10))
    axes = axes.flatten()  # Flatten axes array for easy iteration

    # Columns related to probabilities, for setting y-axis limits later
    probability_columns = ['class_pred_0', 'class_pred_1', 'avg_prob_DNA', 'avg_prob_XNA']

    # Calculate the reference entropy value
    reference_entropy = calculate_shannon_entropy_for_error_rate(correct_p, error_p)

    # Iterate over the columns to be plotted
    for ax, column in zip(axes, plot_columns):
        if column in ['avg_prob_DNA', 'avg_prob_XNA']:
            # Plot with error bars (std deviation) for avg_prob_DNA and avg_prob_XNA
            std_column = f'std_{column[4:]}'  # std_prob_DNA or std_prob_XNA
            summary[column].plot(kind='line', marker='o', ax=ax, yerr=summary[std_column], capsize=5)
        elif column == 'avg_log10_odds_ratio':
            # Plot log10_odds_ratio and add a horizontal line at y=0
            summary[column].plot(kind='line', marker='o', ax=ax, yerr=summary['std_log10_odds_ratio'])
            ax.axhline(y=0, color='r', linestyle='--')
        elif column == 'avg_odds_ratio':
            summary[column].plot(kind='line', marker='o', ax=ax, yerr=summary['std_odds_ratio'])
        elif column == 'avg_first_odds_ratio':
            summary[column].plot(kind='line', marker='o', ax=ax)
        elif column == 'shannon_entropy':
            # Plot Shannon entropy and add the red horizontal line
            summary[column].plot(kind='line', marker='o', ax=ax)
            ax.axhline(y=reference_entropy, color='r', linestyle='--', label=f'Entropy ({correct_p:.2f}/{error_p:.2f})')
            ax.set_ylim(0, 2)  # Ensure the y-axis is set to 2 for entropy
        elif column == 'avg_q_score':
            # Plot avg_q_score with error bars (standard deviation of q_scores)
            summary[column].plot(kind='line', marker='o', ax=ax, yerr=summary['std_q_score'], capsize=5)
        else:
            # Plot the column without error bars
            summary[column].plot(kind='line', marker='o', ax=ax)

        # Map column names to more descriptive titles
        if column == 'class_pred_0':
            title = 'DNA fraction'
        elif column == 'class_pred_1':
            title = 'XNA fraction'
        elif column == 'avg_prob_DNA':
            title = '<P(DNA)>'
        elif column == 'avg_prob_XNA':
            title = '<P(XNA)>'
        elif column == 'avg_log10_odds_ratio':
            title = 'log10(Odds Ratio)'
        elif column == 'avg_odds_ratio':
            title = 'Odds Ratio'
        elif column == 'avg_first_odds_ratio':
            title = 'Average Odds Ratio'
        elif column == 'shannon_entropy':
            title = 'Shannon Entropy'
        elif column == 'avg_q_score':
            title = 'Average Q Score'
        else:
            title = column

        # Set the title and axis labels
        ax.set_title(title)
        ax.set_xlabel('normalized_read_focus_base')
        ax.set_ylabel(title)

        # Set y-axis range from 0 to 1 for probability-related plots
        if column in probability_columns:
            ax.set_ylim(0, 1)

    # Remove any remaining empty subplots if necessary
    for i in range(len(plot_columns), len(axes)):
        fig.delaxes(axes[i])

    # Adjust layout and save the figure
    plt.tight_layout()
    plt.savefig(output_graph_file_path)
    if display_vis:
        plt.show()


def plot_odds_summary(summary, output_odds_graph_file_path):
    """
    Plot the summary statistics using matplotlib.
    
    Parameters:
    - summary: pandas DataFrame, the summary statistics
    
    Returns:
    - None
    """
    # Specify the columns to be plotted
    plot_columns = ['avg_odds_ratio', 'avg_log10_odds_ratio', 'avg_first_odds_ratio'
    ]

     # Set up the subplot grid (1 row, 3 columns) with a larger figure size
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(20, 8))  # Increase figure size
    axes = axes.flatten()  # Flatten axes array for easy iteration

    # Iterate over the columns to be plotted
    for ax, column in zip(axes, plot_columns):
        if column == 'avg_log10_odds_ratio':
            # Plot log10_odds_ratio and add a horizontal line at y=0
            summary[column].plot(kind='line', marker='o', ax=ax, yerr = summary['std_log10_odds_ratio'])
            ax.axhline(y=0, color='r', linestyle='--')
        elif column == 'avg_odds_ratio':
            summary[column].plot(kind='line', marker='o', ax=ax, yerr = summary['std_odds_ratio'] )
        elif column == 'avg_first_odds_ratio':
            summary[column].plot(kind='line', marker='o', ax=ax)
        else:
            # Plot the column without error bars
            summary[column].plot(kind='line', marker='o', ax=ax)

        # Map column names to more descriptive titles
        # Update the title for specific columns
        if column == 'avg_log10_odds_ratio':
            title = 'log10(Odds Ratio)'
        elif column == 'avg_odds_ratio':
            title = 'Odds Ratio'
        elif column == 'avg_first_odds_ratio':
            title = 'Average Odds Ratio'
        else:
            title = column
        
        # Set the title and axis labels
        ax.set_title(title)
        ax.set_xlabel('normalized_read_focus_base')
        ax.set_ylabel(title)

    # Adjust layout and save the figure
    plt.tight_layout()
    plt.savefig(output_odds_graph_file_path)
    if display_vis == True:
        plt.show()

def plot_entropy(summary, save_path, display_vis=False, correct_p=0.85, error_p=0.05):
    """
    Plot the Shannon entropy and base fractions with a horizontal dotted line for a given error rate.
    Stack the base fractions to sum up to 1 for each reference position.
    
    Parameters:
    - summary: pandas DataFrame, the summary statistics
    - save_path: str, path to save the output plot
    - display_vis: bool, whether to display the plot
    - correct_p: float, probability of the correct base being called (e.g., 0.85 for 15% error rate)
    - error_p: float, probability of an incorrect base being called (e.g., 0.05 for 15% error rate)
    
    Returns:
    - None
    """
    # Calculate the reference Shannon entropy for the given error rate
    reference_entropy = calculate_shannon_entropy_for_error_rate(correct_p, error_p)

    # Set up subplots
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(14, 6))
    sns.set_palette("pastel")

    # Plot Shannon entropy
    ax1 = axes[0]
    summary['shannon_entropy'].plot(kind='line', marker='o', ax=ax1, color='b')
    ax1.set_title('Shannon Entropy')
    ax1.set_xlabel('normalized_read_focus_base')
    ax1.set_ylabel('Shannon Entropy')
    ax1.set_ylim(0, 2)  # Explicitly set y-axis limit to 2

    # Adjust x-ticks and labels for the entropy plot
    tick_spacing = 10
    ax1.set_xticks(range(0, len(summary), tick_spacing))
    ax1.set_xticklabels(range(0, len(summary), tick_spacing))

    # Add a dotted horizontal line for the reference Shannon entropy
    ax1.axhline(y=reference_entropy, color='r', linestyle='--', label=f'Theoretical Max Basecalling Entropy ({correct_p:.2f}/{error_p:.2f}*3)')
    ax1.legend()

    # Plot stacked base fractions
    ax2 = axes[1]
    base_columns = ['A', 'C', 'G', 'T']
    
    # Ensure base fractions sum to 1 for each position
    summary_base_fractions = summary[base_columns].copy()
    summary_base_fractions['normalized_read_focus_base'] = summary.index

    # Stacked bar plot: plot each base fraction on top of the others
    summary_base_fractions.set_index('normalized_read_focus_base')[base_columns].plot(
        kind='bar', stacked=True, ax=ax2, color=sns.color_palette("muted"))

    # Set plot labels and title
    ax2.set_title('Stacked Base Fractions')
    ax2.set_xlabel('normalized_read_focus_base')
    ax2.set_ylabel('Fraction')
    ax2.set_ylim(0, 1)  # Ensure the bars stack up to 1
    ax2.set_xticks(range(0, len(summary), tick_spacing))
    ax2.set_xticklabels(range(0, len(summary), tick_spacing))
    
    # Move the legend outside of the plot
    ax2.legend(title="Base", bbox_to_anchor=(1.05, 1), loc='upper left')

    # Add a super title for the entire figure
    plt.suptitle('Shannon Entropy and Stacked Base Fractions', fontsize=16)

    # Adjust layout and save the figure
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(save_path)

    # Display the plot if required
    if display_vis:
        plt.show()



def main():
    """
    Main function to orchestrate the data processing, calculation, and plotting.
    """
    input_file_path = sys.argv[1]
    output_file_path = sys.argv[2]
    graph_file_path = sys.argv[3]
    odds_graph_path = sys.argv[4]
    entropy_graph_file_path = sys.argv[5]

    prob_filter = 0

    df = load_data(input_file_path)
    df_filtered = filter_data(df, prob_filter)
    summary = calculate_summary_statistics(df_filtered)

    save_summary_to_csv(summary, output_file_path)

    plot_summary(summary, graph_file_path)
    plot_odds_summary(summary, odds_graph_path)
    plot_entropy(summary, entropy_graph_file_path)

if __name__ == "__main__":
    main()

