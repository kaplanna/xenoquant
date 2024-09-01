import pandas as pd
import os

def load_csv_file(file_path):
    return pd.read_csv(file_path)

def merge_on_read_id(df1, df2, source_label):
    merged_df = pd.merge(df1, df2, on="read_id", how="inner")
    merged_df["source"] = source_label
    return merged_df

def filter_unmerged_reads(merged_df, all_ids_df):
    merged_read_ids = merged_df["read_id"].unique()
    unmerged_df = all_ids_df[~all_ids_df["read_id"].isin(merged_read_ids)]
    return unmerged_df

def save_csv(df, file_path):
    df.to_csv(file_path, index=False)

# Paths to the input files
all_read_ids_file = "/Users/nickkaplan/xenobiocode/DataAnalysis/basecall/240830_B16_B17_Basecall/demux/all_read_ids.csv"
BA_results_file = "/Users/nickkaplan/xenobiocode/DataAnalysis/basecall/240830_B16_B17_Basecall/BA-Basecall/Results/full_alignment_results.csv"
ST_results_file = "/Users/nickkaplan/xenobiocode/DataAnalysis/basecall/240830_B16_B17_Basecall/ST-Basecall/Results/full_alignment_results.csv"

# Load the CSV files
print("Loading CSV files...")
all_read_ids_df = load_csv_file(all_read_ids_file)
BA_results_df = load_csv_file(BA_results_file)
ST_results_df = load_csv_file(ST_results_file)
print("CSV files loaded.")

# Merge with BA results
print("Merging with BA results...")
merged_df_ba = merge_on_read_id(all_read_ids_df, BA_results_df, "BA")
print(f"BA merge completed. Number of merged reads: {len(merged_df_ba)}")

# Merge with ST results for remaining reads
print("Merging with ST results for remaining reads...")
remaining_reads_df = filter_unmerged_reads(merged_df_ba, all_read_ids_df)
merged_df_st = merge_on_read_id(remaining_reads_df, ST_results_df, "ST")
merged_df = pd.concat([merged_df_ba, merged_df_st])
print(f"ST merge completed. Total number of merged reads: {len(merged_df_st)}")

# Filter out unmerged reads
print("Filtering unmerged reads...")
unmerged_reads_df = filter_unmerged_reads(merged_df, all_read_ids_df)
print(f"Number of unmerged reads: {len(unmerged_reads_df)}")

# Save the merged data and unmerged reads to separate CSV files
output_merged_file = "/Users/nickkaplan/xenobiocode/DataAnalysis/basecall/240830_B16_B17_Basecall/demux/output_merged_reads.csv"
output_unmerged_file = "/Users/nickkaplan/xenobiocode/DataAnalysis/basecall/240830_B16_B17_Basecall/demux/output_unmerged_reads.csv"

print("Saving merged and unmerged reads to CSV files...")
save_csv(merged_df, output_merged_file)
save_csv(unmerged_reads_df, output_unmerged_file)
print("CSV files saved.")

print("Analysis completed.")
