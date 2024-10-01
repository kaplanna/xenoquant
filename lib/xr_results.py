import pandas as pd
import numpy as np
import os
import sys
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score
import xr_params

mod_base = xr_params.mod_base
can_base = xr_params.can_base

# Define directories and paths
def define_directories(base_path):
    results_dir = os.path.join(base_path, 'results')
    bed_dir = os.path.join(base_path, 'references')

    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    
    return results_dir, bed_dir

# Read and process BED file
def read_bed_file(bed_dir, mod_base):
    bed_file_path = os.path.join(bed_dir, f'{mod_base}.bed')
    
    if not os.path.exists(bed_file_path):
        print("Bed file not found in the directory.")
        sys.exit(1)
    
    df_bed = pd.read_csv(bed_file_path, delimiter=r'\s+', header=None)
    df_bed.columns = ['Alignment', 'Start', 'End', 'Name', 'Score', 'Strand']
    return df_bed

# Read the per-read modifications file
def read_data_files(base_path):
    per_read_modifications_file = os.path.join(base_path, 'remora_outputs/per-read_modifications.tsv')
    
    df_modifications = pd.read_csv(per_read_modifications_file, delimiter='\t', dtype={1: 'str', 2: 'str', 3: 'str'})
    df_modifications.columns = ['read_id', 'read_focus_base', 'label', 'class_pred', 'class_probs', 'reference_sequence',
                                'flag', 'ref_start_pos', 'cigar_string', 'ref_length', 'basecalled_sequence', 'q_score']
    return df_modifications

# Process the data (now using reference_sequence instead of sequencing_summary)
def process_data(df_modifications, alignment):
    filtered_rows = df_modifications[df_modifications['reference_sequence'] == alignment]
    
    # Extract class probabilities
    filtered_rows['class_0_probs'] = filtered_rows['class_probs'].str.split(',').str[0].astype(float)
    filtered_rows['class_1_probs'] = 1 - filtered_rows['class_0_probs']
    
    # Ensure class_pred is cast as integer
    filtered_rows['class_pred'] = pd.to_numeric(filtered_rows['class_pred'], errors='coerce').fillna(0).astype(int)
    
    # Retain required columns
    output_df = filtered_rows[['read_id', 'class_pred', 'class_0_probs', 'class_1_probs']].drop_duplicates(subset='read_id')
    return output_df

# Calculate counts and percentages for each alignment
def calculate_results(output_df, alignment):
    number_of_1s = output_df['class_pred'].sum()  # Count 1s (XNA)
    number_of_0s = len(output_df) - number_of_1s  # Count 0s (DNA)
    percentage_1 = round((number_of_1s / len(output_df)) * 100, 2)
    percentage_0 = round((number_of_0s / len(output_df)) * 100, 2)
    
    return {
        'Sequence': alignment,
        'Total Alignments': len(output_df),
        'Number of 1s': number_of_1s,
        'Number of 0s': number_of_0s,
        'Percentage 1': percentage_1,
        'Percentage 0': percentage_0
    }

# Save results to CSV
def save_results(output_df, results_dir, alignment):
    alignment_file_name = f'alignment_results_{alignment}.csv'
    output_file_path = os.path.join(results_dir, alignment_file_name)
    output_df.to_csv(output_file_path, index=False)
    return alignment_file_name

# Main execution function to iterate through alignments and process data
def analyze_results(base_path, plot_confusion_matrix=False, plot_roc_curve=False):
    results_dir, bed_dir = define_directories(base_path)
    df_bed = read_bed_file(bed_dir, mod_base)
    df_modifications = read_data_files(base_path)
    
    full_alignment_results_path = os.path.join(results_dir, 'full_alignment_results.csv')
    if os.path.exists(full_alignment_results_path):
        os.remove(full_alignment_results_path)
    
    results = []
    alignment_files_list = []

    for _, row in df_bed.iterrows():
        alignment = row['Alignment']
        output_df = process_data(df_modifications, alignment)
        output_df.to_csv(full_alignment_results_path, mode='a', header=not os.path.exists(full_alignment_results_path), index=False)
        
        # Calculate and store results
        result = calculate_results(output_df, alignment)
        results.append(result)

        alignment_file_name = save_results(output_df, results_dir, alignment)
        alignment_files_list.append(alignment_file_name)
    
    # Save overall results
    results_df = pd.DataFrame(results)
    final_results_path = os.path.join(results_dir, 'overall_results.csv')
    results_df.to_csv(final_results_path, index=False)
    
    print(f"Full alignment results saved to {full_alignment_results_path}")

    # Optional: plot confusion matrix or ROC curve
    if plot_confusion_matrix:
        plot_conf_matrix(results_df, results_dir)

    if plot_roc_curve:
        plot_roc(alignment_files_list, results_dir)

# Function to plot confusion matrix
def plot_conf_matrix(results_df, results_dir):
    try:
        TP = TN = FP = FN = 0
        xna_type = xr_params.mod_base
        dna_type = xr_params.can_base

        for _, row in results_df.iterrows():
            alignment = row['Sequence'][0]
            if alignment == dna_type:
                FN += row['Number of 1s']
                TN += row['Number of 0s']
            elif alignment == xna_type:
                TP += row['Number of 1s']
                FP += row['Number of 0s']

        conf_matrix_counts = np.array([[TP, FP], [FN, TN]])
        labels = [f"{mod_base}", f"{can_base}"]
        sns.heatmap(conf_matrix_counts, annot=True, fmt='g', cmap='rocket_r', xticklabels=labels, yticklabels=labels)
        plt.xlabel('Predicted')
        plt.ylabel('Actual')
        plt.title('Confusion Matrix (Counts)')
        plt.savefig(os.path.join(results_dir, 'confusion_matrix_counts.pdf'))
        plt.show()
    except Exception as e:
        print(f"Error generating confusion matrix: {e}")

# Function to plot ROC curve
def plot_roc(alignment_files_list, results_dir):
    try:
        y_true_agg = []
        y_prob_agg = []

        for alignment_file in alignment_files_list:
            df_other = pd.read_csv(os.path.join(results_dir, alignment_file))
            true_class = 1 if 'BS' in alignment_file else 0
            y_true_agg.extend([true_class] * len(df_other))
            y_prob_agg.extend(df_other['class_1_probs'])

        fpr, tpr, _ = roc_curve(y_true_agg, y_prob_agg)
        auc = roc_auc_score(y_true_agg, y_prob_agg)
        plt.plot(fpr, tpr, label=f'ROC Curve (AUC = {auc:.2f})', color='orange')
        plt.plot([0, 1], [0, 1], linestyle='--', color='black')
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('ROC Curve')
        plt.legend()
        plt.savefig(os.path.join(results_dir, 'roc_curve.pdf'))
        plt.show()
    except Exception as e:
        print(f"Error generating ROC curve: {e}")

# Run the analysis
if __name__ == "__main__":
    base_path = os.path.expanduser(sys.argv[1])
    analyze_results(base_path, plot_confusion_matrix=False, plot_roc_curve=False)

