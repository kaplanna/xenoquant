########################################################################
########################################################################
"""
xr_results.py 

Title: Unpublished work

By: H. Kawabe, N. Kaplan, J. A. Marchand

Updated: 9/15/23
"""
########################################################################
########################################################################

import pandas as pd
import numpy as np
import os
import sys
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score


# Access the XNA and DNA variables in the params file
import xr_params
mod_base = xr_params.mod_base
can_base = xr_params.can_base


#If you are testing a model where TP/TN/FP/FN is relevant. 
#P / G / A MUST BE THE FIRST LETTERS IN REF FASTA ALIGNMENT NAMES
plot_confusion_matrix = True
plot_roc_curve = True


# Define directories and files based on command line input
base_path = os.path.expanduser(sys.argv[1])
# Check if the base directory exists
if not os.path.exists(base_path):
    print("Xemora [ERROR] - Basecall directory not found. Exiting.")
    sys.exit(1)

# Print status message
print("Xemora [STATUS] - Analyzing Results and Building Plots")

results_dir = os.path.join(base_path, 'Results')
bed_dir = os.path.join(base_path, 'references')

# Create 'Results' directory if it doesn't exist
if not os.path.exists(results_dir):
    os.makedirs(results_dir)

results_dir = os.path.join(base_path, 'Results')
bed_dir = os.path.join(base_path, 'references')

# Create 'Results' directory if it doesn't exist
if not os.path.exists(results_dir):
    os.makedirs(results_dir)

# Define file paths
sequencing_summary_file_txt = os.path.join(base_path, 'preprocess/fastq/sequencing_summary.txt')
per_read_modifications_file = os.path.join(base_path, 'per-read_modifications.tsv')


# Check if either Z.bed or P.bed exists and set the file path accordingly
if os.path.exists(os.path.join(bed_dir, f'{mod_base}.bed')):
    bed_file_path = os.path.join(bed_dir, f'{mod_base}.bed')

else:
    print("Bed file not found in the directory.")
    sys.exit(1)

# Read the .bed file into a Pandas DataFrame
df_bed = pd.read_csv(bed_file_path, delimiter=r'\s+', header=None)
# Name the columns for better readability
df_bed.columns = ['Alignment', 'Start', 'End', 'Name', 'Score', 'Strand']


# Convert sequencing_summary.txt to sequencing_summary.csv
df_sequencing_summary = pd.read_csv(sequencing_summary_file_txt, delimiter='\t')
df_sequencing_summary.to_csv(os.path.join(base_path, 'sequencing_summary.csv'), index=False)

# Read modifications file
df_modifications = pd.read_csv(per_read_modifications_file, delimiter='\t', header=None, dtype={1: 'str', 2: 'str', 3: 'str'})
df_modifications.columns = ['read_id', 'read_focus_base', 'label', 'class_pred', 'class_probs']

# Initialize results dataframe and alignments file names list
results = []
alignment_files_list = []

# Iterate through each row in the .bed file DataFrame
for index, row in df_bed.iterrows():
    sequence_name = row['Name']
    alignment = row['Alignment']
    alignments_file_name = f'alignment_results_{alignment}.csv'
    # Store the alignments_file_name to the list for later use in plotting ROC curves
    alignment_files_list.append(alignments_file_name)

    
    # Filter rows that match the specific alignment
    filtered_rows = df_sequencing_summary[df_sequencing_summary['alignment_genome'] == alignment]
    
    # Create DataFrame for alignment results if it doesn't already exist
    if not os.path.isfile(alignments_file_name):
        df_alignments = pd.DataFrame(columns=['read_id', 'class_pred', 'class_0_probs', 'class_1_probs'])

    
    # Merge based on read_id
    merged_df = pd.merge(filtered_rows, df_modifications, on='read_id')
    
    # Extract class probabilities from the class_probs string
    merged_df['class_0_probs'] = merged_df['class_probs'].str.split(',').str[0].astype(float)
    merged_df['class_1_probs'] = 1 - merged_df['class_0_probs']
    
    # Retain only the columns needed and remove duplicate read_ids
    output_df = merged_df[['read_id', 'class_pred', 'class_0_probs', 'class_1_probs']].drop_duplicates(subset='read_id')
    output_df['class_pred'] = output_df['class_pred'].astype(int)
    
    # Save the filtered and processed data to a CSV file
    output_file_path = os.path.join(results_dir, alignments_file_name)
    output_df.to_csv(output_file_path, index=False)
    
    # Calculate the number of class predictions (1s and 0s)
    number_of_1s = output_df['class_pred'].sum()
    number_of_0s = len(output_df) - number_of_1s
    # Calculate the percentages of 1s and 0s
    percentage_1 = round((number_of_1s / len(output_df)) * 100, 2)
    percentage_0 = round((number_of_0s / len(output_df)) * 100, 2)
    
    # Append calculated data to the results list as a dictionary
    results.append({
        'Sequence': alignment,
        'Total Alignments': len(output_df),
        'Number of 1s': number_of_1s,
        'Number of 0s': number_of_0s,
        'Percentage 1': percentage_1,
        'Percentage 0': percentage_0
    })
# Uncomment to debug and print the 'results' list
# print(results)

# Convert the results list (of dictionaries) to a Pandas DataFrame for easier manipulation and output
results_df = pd.DataFrame(results)

# Save the overall results to a CSV file
print("Xemora [STATUS] - Writing Overall Results File")
final_results_path = os.path.join(results_dir, 'overall_results.csv')
results_df.to_csv(final_results_path, index=False)

# Mapping dictionary
alignment_map = {
    'C': 'G',
    'T': 'A',
    'G': 'G',
    'A': 'A',
    'P': 'P',
    'Z': 'P'
        }
xna_type = alignment_map.get(mod_base, mod_base)
dna_type = alignment_map.get(can_base, can_base)

if plot_confusion_matrix:
    try:
        TP = 0
        FP = 0
        TN = 0
        FN = 0



        for index, row in results_df.iterrows():
            alignment = row['Sequence'][0]

            if alignment == dna_type:
                FN += row['Number of 1s']
                TN += row['Number of 0s']
            elif alignment == xna_type:
                TP += row['Number of 1s']
                FP += row['Number of 0s']

        # Count-based Confusion Matrix
        conf_matrix_counts = np.array([[TP, FP], [FN, TN]])
        labels = [f"{mod_base}", f"{can_base}"]

        sns.set(font_scale=1.2)
        plt.figure(figsize=(8, 6))
        sns.heatmap(conf_matrix_counts, annot=True, fmt='g', cmap='rocket_r', xticklabels=labels, yticklabels=labels)
        plt.xlabel('Predicted')
        plt.ylabel('Actual')
        plt.title('Confusion Matrix (Counts)')
        plt.savefig(os.path.join(results_dir, 'confusion_matrix_counts.pdf'))
        plt.show()

        # Normalized Confusion Matrix
        row_sums = conf_matrix_counts.sum(axis=1, keepdims=True)
        conf_matrix_norm = conf_matrix_counts / row_sums

        plt.figure(figsize=(8, 6))
        sns.heatmap(conf_matrix_norm, annot=True, fmt='g', cmap='rocket_r', xticklabels=labels, yticklabels=labels)
        plt.xlabel('Predicted')
        plt.ylabel('Actual')
        plt.title('Confusion Matrix (Normalized)')
        plt.savefig(os.path.join(results_dir, 'confusion_matrix_normalized.pdf'))
        plt.show()
    except Exception as e:
        print(f"An error occurred: {e}")

if plot_roc_curve == True:


    # Initialize lists to aggregate true class labels and predicted probabilities
    y_true_agg = []
    y_prob_agg = []

    # ... (previous code for reading data and creating alignment-specific CSVs remains unchanged)

    # Iterate through each alignment file to aggregate the data for ROC curve
    for alignments_file_name in alignment_files_list:
        output_file_path = os.path.join(results_dir, alignments_file_name)
        df_other = pd.read_csv(output_file_path)

        # Set the true class based on the alignment
        alignment = alignments_file_name.split('_')[2].split('.')[0]  # Extracting alignment name from file name
        if alignment == xna_type:
            true_class = 1
        elif alignment == dna_type:
            true_class = 0
        else:
            print(f"Unknown alignment: {alignment}")
            continue

        # Append the true class labels and predicted probabilities
        y_true_agg.extend([true_class] * len(df_other))
        y_prob_agg.extend(df_other['class_1_probs'].tolist())

    # Plot the aggregated ROC curve
    try:
        # Compute the ROC curve
        fpr, tpr, thresholds = roc_curve(y_true_agg, y_prob_agg)
        auc = roc_auc_score(y_true_agg, y_prob_agg)
        print(f'The area under ROC curve for {mod_base} vs {can_base} is {auc:.2f}')

        # Plot the ROC curve
        plt.figure(figsize=(10, 8))
        # Draw the ROC curve with specified orange color and line style
        plt.plot(fpr, tpr, label=f'ROC curve for {mod_base} (AUC = {auc:.2f})', linewidth=2, linestyle='-', color='orange')

        # Add a diagonal dashed line with black color
        plt.plot([0, 1], [0, 1], linestyle='--', color='black')

        # Customize axis scale
        plt.axis('square')  # Make the plot square
        plt.xlim([0, 1])
        plt.ylim([0, 1])
        plt.gca().add_patch(plt.Rectangle((0, 0), 1, 1, fill=False, edgecolor='black', linewidth=2))
        plt.gca().set_facecolor('white')
        plt.xlabel('False Positive Rate', fontsize=14)
        plt.ylabel('True Positive Rate', fontsize=14)
        plt.title('Receiver Operating Characteristic (ROC) Curve', fontsize=16)
        plt.grid(False)
        plt.legend(loc='lower right', fontsize=12)
        plt.savefig(os.path.join(results_dir, 'roc_curve.pdf'))
        plt.show()

    except Exception as e:
        print(f"Failed to generate ROC curve. Error: {e}")
        
 





