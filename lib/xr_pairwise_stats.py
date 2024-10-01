########################################################################
########################################################################
"""
xr_pairwise_stats.py 

Title: Unpublished work

Calculates pairwise analysis from 2 tsv files containing consensus data

[Logs]
9/18 - file generated

By: H. Kawabe, N. Kaplan, J. Sumabat, J. A. Marchand

Updated: 11/28/23
"""
########################################################################
########################################################################

import os
import glob
import sys
import tempfile
import matplotlib.pyplot as plt
import seaborn as sns

from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path
from xr_tools  import *
from xr_params import *
from collections import defaultdict

########################################################################
########################################################################
output_dir = sys.argv[1]
fwd_data = sys.argv[2]
rev_data = sys.argv[3]

########################################################################
########################################################################

#Initialize Directories
output_dir = check_make_dir(output_dir)

# Function to load data from forward and reverse files
def load_data(fwd_file, rev_file):
    """
    Loads two dataframes from forward and reverse files and merges them horizontally
    based on the 'normalized_read_focus_base' column, after reversing the base order
    for the reverse file.

    Parameters:
    - fwd_file: str, path to the forward data file (TSV format)
    - rev_file: str, path to the reverse data file (TSV format)

    Returns:
    - merged_df: pandas DataFrame, the horizontally merged dataframe with suffixes
    """
    
    print('Xemora [STATUS] - Loading and merging consensus results for pairwise analysis')

    # Load forward and reverse DataFrames
    fwd_df = pd.read_csv(fwd_file, sep='\t')
    rev_df = pd.read_csv(rev_file, sep='\t')

    # Find the maximum base value from both dataframes
    max_base = max(fwd_df['normalized_read_focus_base'].max(), rev_df['normalized_read_focus_base'].max())
    min_base = fwd_df['normalized_read_focus_base'].min()  # Get minimum base value

    # Reverse the 'normalized_read_focus_base' column for reverse dataframe
    rev_df['normalized_read_focus_base'] = max_base - rev_df['normalized_read_focus_base'] + min_base

    # Merge the DataFrames horizontally on 'normalized_read_focus_base' with suffixes
    merged_df = pd.merge(fwd_df, rev_df, on='normalized_read_focus_base', suffixes=('_fwd', '_rev'))

    return merged_df

def calculate_pairwise_stats(merged_df):
    print('Xemora [STATUS] - Calculating pairwise stats')

    # Calculating average XNA and DNA recall
    merged_df['class_pred_0_pairwise'] = (merged_df['class_pred_0_fwd'] + merged_df['class_pred_0_rev'])/2
    merged_df['class_pred_1_pairwise'] = (merged_df['class_pred_1_fwd'] + merged_df['class_pred_1_rev'])/2

    # Calculate geometric mean for DNA and XNA probabilities (2 strand system)
    merged_df['avg_prob_XNA_pairwise_geo_mean'] = np.sqrt(merged_df['avg_prob_XNA_fwd'] * merged_df['avg_prob_XNA_rev'])
    merged_df['avg_prob_DNA_pairwise_geo_mean'] = np.sqrt(merged_df['avg_prob_DNA_fwd'] * merged_df['avg_prob_DNA_rev'])

    # Calculate error from XNA and DNA geometric means (2 strand system)
    merged_df['avg_prob_XNA_pairwise_error'] = (merged_df['avg_prob_XNA_pairwise_geo_mean'] / 2 * np.sqrt((merged_df['std_prob_XNA_fwd'] / merged_df['avg_prob_XNA_fwd']) ** 2 +(merged_df['std_prob_XNA_rev'] / merged_df['avg_prob_XNA_rev']) ** 2))
    merged_df['avg_prob_DNA_pairwise_error'] = (merged_df['avg_prob_DNA_pairwise_geo_mean'] / 2 * np.sqrt((merged_df['std_prob_DNA_fwd'] / merged_df['avg_prob_DNA_fwd']) ** 2 +(merged_df['std_prob_DNA_rev'] / merged_df['avg_prob_DNA_rev']) ** 2))

    # Calculate the sum of the forward and reverse avg_log10_odds_ratio columns
    '''
    might need to divide by 2 here
    '''
    merged_df['avg_log10_odds_ratio_pairwise_sum'] = (merged_df['avg_log10_odds_ratio_fwd'] + merged_df['avg_log10_odds_ratio_rev'])
    
    # Calculate error from avg_log10_odds_ratio standard deviations
    merged_df['avg_log10_odds_ratio_error'] = np.sqrt(merged_df['std_log10_odds_ratio_fwd']**2 + merged_df['std_log10_odds_ratio_rev']**2)

    #  Shannon entropy here maybe 

    # Calculate consensus quality score
    merged_df['avg_q_score_pairwise'] = (merged_df['avg_q_score_fwd'] + merged_df['avg_q_score_rev'])/2

    # Calculate error for pairwise q score
    merged_df['std_q_score_pairwise'] = (np.sqrt((merged_df['std_q_score_fwd'] / 2) ** 2 +(merged_df['std_q_score_rev'] / 2) ** 2))
    return merged_df

def plot_pairwise_stats(merged_df, output_dir, overall_title="Pairwise and Strand Statistics", display_vis=False):
    """
    Plots pairwise statistics along with individual forward and reverse strand statistics
    as a function of 'normalized_read_focus_base', with consistent colors, better visibility, 
    a single shared legend, customizable subplot titles and y-axis labels, and an overall title.

    Parameters:
    - merged_df: pandas DataFrame containing the merged forward and reverse data.
    - output_dir: str, path to the directory where the plots will be saved.
    - overall_title: str, the overall title for the entire figure (shown above the legend).
    - display_vis: bool, whether to display the plot interactively (default is False).
    """
    
    # Set seaborn style directly and use a colorblind-friendly palette
    sns.set(style="darkgrid", context="paper")
    palette = sns.color_palette("Set2", 3)  # Softer color palette for forward, reverse, and pairwise

    # Custom titles and y-axis labels for each plot
    titles_and_labels = {
        'class_pred_0_fwd': ('DNA Class Prediction', 'Fraction of DNA Reads'),
        'class_pred_1_fwd': ('XNA Class Prediction', 'Fraction of XNA Reads'),
        'avg_prob_XNA_fwd': ('Average XNA Probability', 'Probability'),
        'avg_prob_DNA_fwd': ('Average DNA Probability', 'Probability'),
        'avg_log10_odds_ratio_fwd': ('Log10 Odds Ratio', 'Log10(Odds Ratio)'),
        'avg_q_score_fwd': ('Average Q-Score', 'Q-Score'),
        'avg_first_odds_ratio_fwd': ('First Odds Ratio', 'Odds Ratio'),
        'shannon_entropy_fwd': ('Shannon Entropy', 'Entropy'),
        'avg_odds_ratio_fwd': ('Average Odds Ratio', 'Odds Ratio')
    }

    # Define the columns to be plotted for forward, reverse, and pairwise data
    pairwise_columns_to_plot = {
        'class_pred_0_fwd': ('class_pred_0_rev', 'class_pred_0_pairwise'),
        'class_pred_1_fwd': ('class_pred_1_rev', 'class_pred_1_pairwise'),
        'avg_prob_XNA_fwd': ('avg_prob_XNA_rev', 'avg_prob_XNA_pairwise_geo_mean'),
        'avg_prob_DNA_fwd': ('avg_prob_DNA_rev', 'avg_prob_DNA_pairwise_geo_mean'),
        'avg_log10_odds_ratio_fwd': ('avg_log10_odds_ratio_rev', 'avg_log10_odds_ratio_pairwise_sum'),
        'avg_q_score_fwd': ('avg_q_score_rev', 'avg_q_score_pairwise')
    }

    # Error columns for forward, reverse, and pairwise data
    error_columns = {
        'avg_prob_XNA_fwd': ('std_prob_XNA_fwd', 'std_prob_XNA_rev', 'avg_prob_XNA_pairwise_error'),
        'avg_prob_DNA_fwd': ('std_prob_DNA_fwd', 'std_prob_DNA_rev', 'avg_prob_DNA_pairwise_error'),
        'avg_log10_odds_ratio_fwd': ('std_log10_odds_ratio_fwd', 'std_log10_odds_ratio_rev', 'avg_log10_odds_ratio_error'),
        'avg_q_score_fwd': ('std_q_score_fwd', 'std_q_score_rev', 'std_q_score_pairwise')
    }

    # Additional columns that don't have pairwise analysis
    additional_columns_to_plot = [
        'avg_first_odds_ratio_fwd',  # For forward and reverse only
        'shannon_entropy_fwd',       # For forward and reverse only
        'avg_odds_ratio_fwd'         # For forward and reverse only
    ]

    # Corresponding reverse columns for the additional plots
    additional_reverse_columns = [
        'avg_first_odds_ratio_rev',
        'shannon_entropy_rev',
        'avg_odds_ratio_rev'
    ]

    # Create a figure with a 3x3 layout (9 total plots)
    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(18, 18))
    axes = axes.flatten()  # Flatten to easily iterate

    # Colors: forward, reverse, and pairwise
    forward_color = palette[0]
    reverse_color = palette[1]
    pairwise_color = palette[2]

    # Store handles and labels for legend
    handles = []
    labels = []

    # Plot each pairwise, forward, and reverse column
    for i, (fwd_col, (rev_col, pairwise_col)) in enumerate(pairwise_columns_to_plot.items()):
        ax = axes[i]
        
        # Get the custom title and y-label for the plot
        plot_title, y_label = titles_and_labels.get(fwd_col, (fwd_col, fwd_col))
        
        ax.set_title(plot_title, fontsize=12)
        ax.set_xlabel('normalized_read_focus_base', fontsize=10)
        ax.set_ylabel(y_label, fontsize=10)
        
        # Plot the forward strand data
        line_fwd, = ax.plot(merged_df['normalized_read_focus_base'], merged_df[fwd_col], 'o-', color=forward_color, label='Forward Strand')

        # Plot the reverse strand data
        line_rev, = ax.plot(merged_df['normalized_read_focus_base'], merged_df[rev_col], 'x-', color=reverse_color, label='Reverse Strand')

        # Plot the pairwise data
        line_pairwise, = ax.plot(merged_df['normalized_read_focus_base'], merged_df[pairwise_col], 's-', color=pairwise_color, label='Pairwise Analysis')

        # Add error bars if applicable
        if fwd_col in error_columns:
            fwd_error, rev_error, pairwise_error = error_columns[fwd_col]
            ax.errorbar(merged_df['normalized_read_focus_base'], merged_df[fwd_col], 
                        yerr=merged_df[fwd_error], fmt='o-', color=forward_color, capsize=5)
            ax.errorbar(merged_df['normalized_read_focus_base'], merged_df[rev_col], 
                        yerr=merged_df[rev_error], fmt='x-', color=reverse_color, capsize=5)
            ax.errorbar(merged_df['normalized_read_focus_base'], merged_df[pairwise_col], 
                        yerr=merged_df[pairwise_error], fmt='s-', color=pairwise_color, capsize=5)

        # Set y-axis limit for probabilities and class predictions (max at 1)
        if 'class_pred' in fwd_col or 'avg_prob' in fwd_col:
            ax.set_ylim(0, 1)

        # Ensure no overlapping labels
        ax.tick_params(axis='x', rotation=45)

        # Add handles and labels for the legend (only once)
        if i == 0:
            handles.extend([line_fwd, line_rev, line_pairwise])
            labels.extend(['Forward Strand', 'Reverse Strand', 'Pairwise Analysis'])

    # Plot the remaining original (non-pairwise) forward and reverse plots
    for j, (fwd_col, rev_col) in enumerate(zip(additional_columns_to_plot, additional_reverse_columns), start=i + 1):
        ax = axes[j]

        # Get the custom title and y-label for the plot
        plot_title, y_label = titles_and_labels.get(fwd_col, (fwd_col, fwd_col))
        
        ax.set_title(plot_title, fontsize=12)
        ax.set_xlabel('normalized_read_focus_base', fontsize=10)
        ax.set_ylabel(y_label, fontsize=10)

        # Plot the forward strand data
        ax.plot(merged_df['normalized_read_focus_base'], merged_df[fwd_col], 'o-', color=forward_color, label='Forward Strand')

        # Plot the reverse strand data
        ax.plot(merged_df['normalized_read_focus_base'], merged_df[rev_col], 'x-', color=reverse_color, label='Reverse Strand')

        # Add error bars if applicable
        if fwd_col == 'avg_odds_ratio_fwd':
            ax.errorbar(merged_df['normalized_read_focus_base'], merged_df[fwd_col], 
                        yerr=merged_df['std_odds_ratio_fwd'], fmt='o-', color=forward_color, capsize=5)
            ax.errorbar(merged_df['normalized_read_focus_base'], merged_df[rev_col], 
                        yerr=merged_df['std_odds_ratio_rev'], fmt='x-', color=reverse_color, capsize=5)

        # Set y-axis limit where appropriate
        if fwd_col in ['avg_prob_DNA_fwd', 'avg_prob_XNA_fwd']:
            ax.set_ylim(0, 1)

    # Adjust layout to ensure everything fits well and give more space for the legend
    plt.tight_layout(rect=[0, 0, 1, 0.93])  # Leave more space at the top for title and legend

    # Add overall title above the legend
    fig.suptitle(overall_title, fontsize=16)

    # Add a single legend just below the overall title
    fig.legend(handles=handles, labels=labels, loc='upper center', ncol=3, bbox_to_anchor=(0.5, 0.97), fontsize=12)

    # Save the figure as PDF
    pdf_output = os.path.join(output_dir, 'forward_reverse_pairwise_stats_plots.pdf')
    plt.savefig(pdf_output, bbox_inches='tight')

    # Optionally show the plot
    if display_vis:
        plt.show()


def main():
    
    merged_df = load_data(fwd_data, rev_data)
    merged_df = calculate_pairwise_stats(merged_df)
    merged_df.to_csv(os.path.join(output_dir, 'merged_df.tsv'), sep='\t', index=False)
    plot_pairwise_stats(merged_df, output_dir)
if __name__ == "__main__":
    main()
