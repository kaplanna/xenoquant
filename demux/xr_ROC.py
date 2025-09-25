import pandas as pd
from sklearn.metrics import roc_curve, roc_auc_score
import matplotlib.pyplot as plt
import os

# Load the BA and ST results TSV files with the correct delimiter
BA_results_df = pd.read_csv(
    '/home/marchandlab/DataAnalysis/Kaplan/basecall/DsPx/250904_D12_Basecall/DsN-Basecall/demux/demux_per-read_modifications.tsv',
    sep='\t'
)
ST_results_df = pd.read_csv(
    '/home/marchandlab/DataAnalysis/Kaplan/basecall/DsPx/250904_D12_Basecall/PxN-Basecall/demux/demux_per-read_modifications.tsv',
    sep='\t'
)
# Specify the output directory for saving the plot
output_dir = '/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/DsPx/250214_D1_Basecall/ROC_Curves'  # Replace with your desired directory
os.makedirs(output_dir, exist_ok=True)  # Ensure the directory exists
# File name for the saved plot
output_file = os.path.join(output_dir, 'DsPx_ROC_Curve.pdf')


# Select only the required columns
BA_results_df = BA_results_df[['read_id', 'class_pred', 'class_probs', 'barcode_pair']]
ST_results_df = ST_results_df[['read_id', 'class_pred', 'class_probs', 'barcode_pair']]
# Extract class_1_prob from class_probs
BA_results_df['class_1_prob'] = BA_results_df['class_probs'].str.split(',').str[1].astype(float)
ST_results_df['class_1_prob'] = ST_results_df['class_probs'].str.split(',').str[1].astype(float)

# Assign ground truth labels
ground_truth_map = {
    'NB12_FWD_NB31_REV': 0,
    'None_FWD_NB30_REV': 1
}

BA_results_df['ground_truth'] = BA_results_df['barcode_pair'].map(ground_truth_map)
ST_results_df['ground_truth'] = ST_results_df['barcode_pair'].map(ground_truth_map)

# Filter out rows where ground_truth is NaN
BA_results_df = BA_results_df.dropna(subset=['ground_truth'])
ST_results_df = ST_results_df.dropna(subset=['ground_truth'])

# Generate ROC curve for BA
fpr_ba, tpr_ba, _ = roc_curve(BA_results_df['ground_truth'], BA_results_df['class_1_prob'])
roc_auc_ba = roc_auc_score(BA_results_df['ground_truth'], BA_results_df['class_1_prob'])

# Generate ROC curve for ST
fpr_st, tpr_st, _ = roc_curve(ST_results_df['ground_truth'], ST_results_df['class_1_prob'])
roc_auc_st = roc_auc_score(ST_results_df['ground_truth'], ST_results_df['class_1_prob'])

# Plot both ROC curves on the same figure
plt.figure(figsize=(3, 3), dpi=300)  # Smaller figure size with high resolution
plt.plot(
    fpr_ba, tpr_ba, color='#375EDB', lw=1, 
    label=f'Ds-N ROC (AUC = {roc_auc_ba:.2f})'
)  # Thin line for BA
plt.plot(
    fpr_st, tpr_st, color='#31DE80', lw=1, 
    label=f'Px-N ROC (AUC = {roc_auc_st:.2f})'
)  # Thin line for ST
plt.plot([0, 1], [0, 1], color='black', lw=0.5, linestyle='--')  # Thinner diagonal line

plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xlabel('False Positive Rate', fontsize=8, fontname='Arial')
plt.ylabel('True Positive Rate', fontsize=8, fontname='Arial')
plt.title('Receiver Operating Characteristic', fontsize=10, fontweight='bold', fontname='Arial')
plt.legend(loc='lower right', fontsize=6, frameon=True, edgecolor='black')
plt.gca().set_aspect('equal', adjustable='box')
plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray')
plt.gca().spines['top'].set_visible(True)
plt.gca().spines['right'].set_visible(True)
plt.gca().spines['left'].set_color('black')
plt.gca().spines['bottom'].set_color('black')
plt.gca().tick_params(width=0.5, labelsize=6)  # Smaller tick size
plt.tight_layout()  # Ensures no clipping

# Save the plot to the specified directory
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"Plot saved to {output_file}")


plt.show()


