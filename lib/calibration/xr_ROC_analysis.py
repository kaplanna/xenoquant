import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, precision_recall_curve, f1_score
import os

### USER SELECTION: Choose which threshold to use ###
use_f1 = False  # Set to False to use Youden’s J instead of F1-max

# Define barcode pairs for true positive and true negative samples
true_positive_barcode = "NB25_FWD_NB25_REV"
true_negative_barcode = "NB14_FWD_NB10_REV"

# Load the calibrated dataset
demux_dir = "/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/PZn/250213_P1_Basecall/ZC_Basecall/demux"
file_path = os.path.join(demux_dir, "recalibrated_PRM.csv")
df = pd.read_csv(file_path)

# Filter dataset to use only the ground truth positive and negative barcode pairs
df_pos = df[df['barcode_pair'] == true_positive_barcode].copy()
df_neg = df[df['barcode_pair'] == true_negative_barcode].copy()

# Assign ground truth labels
df_pos['true_label'] = 1  # True positives
df_neg['true_label'] = 0  # True negatives

# Combine for ROC and F1 analysis
df_analysis = pd.concat([df_pos, df_neg], ignore_index=True)

# Extract true labels and calibrated probabilities
y_true = df_analysis['true_label'].values
y_scores = df_analysis['calibrated_prob'].values

# ------------------------- Compute ROC Curve & Youden’s J Threshold -------------------------
# Compute ROC curve and AUC
fpr, tpr, roc_thresholds = roc_curve(y_true, y_scores)
roc_auc = auc(fpr, tpr)

# Find the best classification threshold using Youden’s J statistic (maximizing TPR - FPR)
youden_index = np.argmax(tpr - fpr)
best_youden_threshold = roc_thresholds[youden_index]

# Plot ROC curve (for reference)
plt.figure(figsize=(8, 6))
plt.plot(fpr, tpr, label=f'ROC Curve (AUC = {roc_auc:.3f})', color='blue')
plt.scatter(fpr[youden_index], tpr[youden_index], color='red', label=f"Best Youden Threshold: {best_youden_threshold:.3f}", zorder=3)
plt.plot([0, 1], [0, 1], 'k--', label="Random Classifier (AUC = 0.5)")

plt.xlabel("False Positive Rate (FPR)")
plt.ylabel("True Positive Rate (TPR)")
plt.title("ROC Curve for Reference")
plt.legend()
plt.grid()
fig_path = os.path.join(demux_dir, "ROC_Cal.pdf")
plt.savefig(fig_path)
plt.show()

# ------------------------- Compute Precision-Recall Curve & Find F1-max Threshold -------------------------
# Compute Precision-Recall curve
precisions, recalls, pr_thresholds = precision_recall_curve(y_true, y_scores)

# Compute F1 scores for each threshold
f1_scores = 2 * (precisions * recalls) / (precisions + recalls + 1e-10)  # Avoid division by zero

# Find the threshold that maximizes the F1-score
best_f1_index = np.argmax(f1_scores)
best_f1_threshold = pr_thresholds[best_f1_index]
best_f1_score = f1_scores[best_f1_index]

# Plot Precision-Recall Curve with F1-max Threshold
plt.figure(figsize=(8, 6))
plt.plot(pr_thresholds, f1_scores[:-1], label="F1-score", color="blue")
plt.axvline(best_f1_threshold, linestyle="--", color="red", label=f"Best F1 Threshold: {best_f1_threshold:.3f}")
plt.xlabel("Probability Threshold")
plt.ylabel("F1-score")
plt.title("F1-score vs. Classification Threshold")
plt.legend()
plt.grid()
f1_fig_path = os.path.join(demux_dir, "F1_Cal.pdf")
plt.savefig(f1_fig_path)
plt.show()

# ------------------------- Choose the Classification Threshold Based on `use_f1` -------------------------
if use_f1:
    selected_threshold = best_f1_threshold
    classification_method = "F1-max"
else:
    selected_threshold = best_youden_threshold
    classification_method = "Youden’s J"

# Apply the chosen threshold to classify all reads
df['final_class'] = (df['calibrated_prob'] >= selected_threshold).astype(int)

# Save the updated dataset with new classifications
output_path = os.path.join(demux_dir, "reclassified_PRM.csv")
df.to_csv(output_path, index=False)

# ------------------------- Print Summary -------------------------
print(f"\nSelected Classification Method: {classification_method}")
print(f"Optimal {classification_method} Threshold: {selected_threshold:.3f}")
print(f"AUC Score for Reference: {roc_auc:.3f}")
print(f"Updated classifications saved in '{output_path}'.")


