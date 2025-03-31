import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression
from sklearn.calibration import calibration_curve
from sklearn.model_selection import train_test_split
import joblib
import os

### User-defined barcode pairs for calibration ###
true_positive_barcode = "NB25_FWD_NB25_REV"  # Barcode where all reads are true positives (ground truth = 1)
true_negative_barcode = "NB14_FWD_NB10_REV"  # Barcode where all reads are true negatives (ground truth = 0)


demux_dir = "/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/PZn/250213_P1_Basecall/PG_Basecall/demux"
file_path = os.path.join(demux_dir, "demux_per-read_modifications.tsv")


df = pd.read_csv(file_path, sep='\t')

# Extract class probabilities (parsing "class_probs" column)
df[['class_0_prob', 'class_1_prob']] = df['class_probs'].str.split(',', expand=True).astype(float)

# Identify true positives and true negatives based on barcode pairs
df_pos = df[df['barcode_pair'] == true_positive_barcode].copy()
df_neg = df[df['barcode_pair'] == true_negative_barcode].copy()
df_pos['true_label'] = 1
df_neg['true_label'] = 0

# Combine true positives and negatives for calibration
df_calibration = pd.concat([df_pos, df_neg], ignore_index=True)

# Extract probabilities and true labels for training calibration models
pred_probs = df_calibration['class_1_prob'].values  # ML confidence scores
true_labels = df_calibration['true_label'].values  # True class labels

# Split into training and test sets for calibration
X_train, X_test, y_train, y_test = train_test_split(pred_probs.reshape(-1, 1), true_labels, test_size=0.2, random_state=42)

# Platt Scaling (Logistic Regression)
platt_model = LogisticRegression()
platt_model.fit(X_train, y_train)
platt_probs = platt_model.predict_proba(X_test)[:, 1]  # Get new probabilities

# Isotonic Regression
iso_model = IsotonicRegression(out_of_bounds='clip')
iso_model.fit(X_train.ravel(), y_train)
iso_probs = iso_model.predict(X_test.ravel())

# Compare calibration curves
plt.figure(figsize=(8, 6))
plt.plot([0, 1], [0, 1], "k--", label="Perfect Calibration (y=x)")

# Before calibration
fraction_of_positives, mean_predicted_value = calibration_curve(y_test, X_test.ravel(), n_bins=10)
plt.plot(mean_predicted_value, fraction_of_positives, "s-", label="Before Calibration")

# Platt Scaling calibration
fraction_of_positives_platt, mean_predicted_value_platt = calibration_curve(y_test, platt_probs, n_bins=10)
plt.plot(mean_predicted_value_platt, fraction_of_positives_platt, "o-", label="Platt Scaling")

# Isotonic Regression calibration
fraction_of_positives_iso, mean_predicted_value_iso = calibration_curve(y_test, iso_probs, n_bins=10)
plt.plot(mean_predicted_value_iso, fraction_of_positives_iso, "d-", label="Isotonic Regression")

plt.xlabel("Mean Predicted Probability")
plt.ylabel("Fraction of True 1s")
plt.title("Calibration Curve (Using True Positive & True Negative Barcodes)")
plt.legend()
plt.grid()
fig_path = os.path.join(demux_dir, "calibration_plot.pdf")
plt.savefig(fig_path)
plt.show()

# Determine the best calibration method (defaulting to Isotonic Regression if it looks best)
best_model = iso_model
print("Using Isotonic Regression for recalibration.")


# Apply the best calibration model to the entire dataset (not just true pos/neg barcodes)
df['calibrated_prob'] = best_model.predict(df['class_1_prob'].values)

# Save the recalibrated dataset
output_path = os.path.join(demux_dir, "recalibrated_PRM.csv")
df.to_csv(output_path, index=False)

print("Recalibration complete! All reads have been adjusted and saved.")

