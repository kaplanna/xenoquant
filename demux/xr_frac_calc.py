import pandas as pd

def load_csv_file(file_path):
    return pd.read_csv(file_path)

def calculate_xna_fractions(merged_df):
    # Separate the BA and ST records
    ba_df = merged_df[merged_df['source'] == 'BA']
    st_df = merged_df[merged_df['source'] == 'ST']

    # Group by barcode_pair and class_pred to count occurrences for BA
    ba_tally_df = ba_df.groupby(['barcode_pair', 'class_pred']).size().unstack(fill_value=0)
    ba_tally_df['Total_Sequences_BA'] = ba_tally_df.sum(axis=1)
    ba_tally_df['Num_1s_BA'] = ba_tally_df.get(1, 0)
    ba_tally_df['Num_0s_BA'] = ba_tally_df.get(0, 0)
    ba_tally_df['XNA_fraction_B'] = ba_tally_df['Num_1s_BA'] / ba_tally_df['Total_Sequences_BA']

    # Group by barcode_pair and class_pred to count occurrences for ST
    st_tally_df = st_df.groupby(['barcode_pair', 'class_pred']).size().unstack(fill_value=0)
    st_tally_df['Total_Sequences_ST'] = st_tally_df.sum(axis=1)
    st_tally_df['Num_1s_ST'] = st_tally_df.get(1, 0)
    st_tally_df['Num_0s_ST'] = st_tally_df.get(0, 0)
    st_tally_df['XNA_fraction_S'] = st_tally_df['Num_1s_ST'] / st_tally_df['Total_Sequences_ST']

    # Combine results into a single DataFrame
    tally_df = pd.merge(
        ba_tally_df[['Total_Sequences_BA', 'Num_1s_BA', 'Num_0s_BA', 'XNA_fraction_B']],
        st_tally_df[['Total_Sequences_ST', 'Num_1s_ST', 'Num_0s_ST', 'XNA_fraction_S']],
        left_index=True, right_index=True, how='outer'
    )

    # Fill any missing values with 0 (this would occur if a barcode pair is present in only one source)
    tally_df = tally_df.fillna(0)

    # Reset index to bring 'barcode_pair' as a column for splitting
    tally_df = tally_df.reset_index()

    # Split the 'barcode_pair' into 'FWD_barcode' and 'REV_barcode'
    barcode_split = tally_df['barcode_pair'].str.split('_', expand=True)
    tally_df['FWD_barcode'] = barcode_split[0]
    tally_df['REV_barcode'] = barcode_split[2]

    # Reorder the columns
    tally_df = tally_df[['FWD_barcode', 'REV_barcode', 'Total_Sequences_BA', 'Num_1s_BA', 'Num_0s_BA', 'XNA_fraction_B', 
                         'Total_Sequences_ST', 'Num_1s_ST', 'Num_0s_ST', 'XNA_fraction_S']]

    return tally_df


def save_csv(df, file_path):
    df.to_csv(file_path, index=False)

# Paths to the input and output files
merged_results_file = "/Users/nickkaplan/DataAnalysis/basecall/240830_B16_B17_Basecall/demux/output_merged_reads.csv"
output_tally_file = "/Users/nickkaplan/DataAnalysis/basecall/240830_B16_B17_Basecall/demux/output_xna_fractions.csv"

# Load the merged results CSV file
print("Loading merged results CSV file...")
merged_df = load_csv_file(merged_results_file)
print("Merged results CSV file loaded.")

# Calculate XNA fractions with detailed stats
print("Calculating XNA fractions with detailed stats...")
xna_fractions_df = calculate_xna_fractions(merged_df)
print("XNA fractions calculated.")

# Save the results to a new CSV file
print("Saving XNA fractions to CSV file...")
save_csv(xna_fractions_df, output_tally_file)
print("XNA fractions CSV file saved.")

print("XNA fraction analysis with details completed.")
