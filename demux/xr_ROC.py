import pandas as pd
from sklearn.metrics import roc_curve, roc_auc_score
import matplotlib.pyplot as plt

# Load the CSV files
all_read_ids_df = pd.read_csv('/Users/nickkaplan/DataAnalysis/basecall/240830_B16_B17_Basecall/demux/all_read_ids.csv')
BA_results_df = pd.read_csv('/Users/nickkaplan/DataAnalysis/basecall/240830_B16_B17_Basecall/BA-Basecall/Results/full_alignment_results.csv')
ST_results_df = pd.read_csv('/Users/nickkaplan/DataAnalysis/basecall/240830_B16_B17_Basecall/ST-Basecall/Results/full_alignment_results.csv')

# Merge the read IDs with BA and ST results
BA_merged_df = pd.merge(all_read_ids_df, BA_results_df, on='read_id', how='inner')
ST_merged_df = pd.merge(all_read_ids_df, ST_results_df, on='read_id', how='inner')

# Assign ground truth labels
ground_truth_map = {
    'NB12_FWD_NB20_REV': 0,
    'NB25_FWD_NB25_REV': 1
}

BA_merged_df['ground_truth'] = BA_merged_df['barcode_pair'].map(ground_truth_map)
ST_merged_df['ground_truth'] = ST_merged_df['barcode_pair'].map(ground_truth_map)

# Filter out rows where ground_truth is NaN (if any)
BA_merged_df = BA_merged_df.dropna(subset=['ground_truth'])
ST_merged_df = ST_merged_df.dropna(subset=['ground_truth'])

# Generate ROC curve for BA
fpr_ba, tpr_ba, _ = roc_curve(BA_merged_df['ground_truth'], BA_merged_df['class_1_probs'])
roc_auc_ba = roc_auc_score(BA_merged_df['ground_truth'], BA_merged_df['class_1_probs'])

# Generate ROC curve for ST
fpr_st, tpr_st, _ = roc_curve(ST_merged_df['ground_truth'], ST_merged_df['class_1_probs'])
roc_auc_st = roc_auc_score(ST_merged_df['ground_truth'], ST_merged_df['class_1_probs'])

# Plot both ROC curves on the same figure
plt.figure(figsize=(6, 6))
plt.plot(fpr_ba, tpr_ba, color='#375EDB', lw=3, label=f'BA ROC (AUC = {roc_auc_ba:.2f})', marker='o', markersize=5)
plt.plot(fpr_st, tpr_st, color='#31DE80', lw=3, label=f'ST ROC (AUC = {roc_auc_st:.2f})', marker='o', markersize=5)
plt.plot([0, 1], [0, 1], color='black', lw=2, linestyle='--')

plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.xlabel('False Positive Rate', fontsize=14, fontname='Arial')
plt.ylabel('True Positive Rate', fontsize=14, fontname='Arial')
plt.title('Receiver Operating Characteristic', fontsize=16, fontweight='bold', fontname='Arial')
plt.legend(loc='lower right', fontsize=12, frameon=True, edgecolor='black')
plt.gca().set_aspect('equal', adjustable='box')
plt.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray')
plt.gca().spines['top'].set_visible(True)
plt.gca().spines['right'].set_visible(True)
plt.gca().spines['left'].set_color('black')
plt.gca().spines['bottom'].set_color('black')
plt.gca().tick_params(width=2)
plt.show()
