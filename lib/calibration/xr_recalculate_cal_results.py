import pandas as pd
import os

# Load the reclassified dataset (now using F1-max classification)
demux_dir = "/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/PZn/250213_P1_Basecall/ZC_Basecall/demux"
file_path = os.path.join(demux_dir, "reclassified_PRM.csv")
df = pd.read_csv(file_path)

# Group by barcode pair and compute summary statistics using `final_class`
summary_df = df.groupby("barcode_pair").agg(
    number_of_alignments=("final_class", "count"),  # Total occurrences per barcode pair
    number_of_1s=("final_class", "sum")  # Count of ones per barcode pair
).reset_index()

# Compute fraction of 1s
summary_df["fraction_1"] = summary_df["number_of_1s"] / summary_df["number_of_alignments"]

# Save the summary file
output_file = os.path.join(demux_dir, "reclassified_SUMMARY.csv")
summary_df.to_csv(output_file, index=False)

print(f"Summary saved as '{output_file}'.")

