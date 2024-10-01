########################################################################
########################################################################
"""
xr_realign.py 

Title: Unpublished work

xr_realign.py performs normalization and realignment and exracts useful columns 
for plotting

To do: add more print statements


By: J. Sumabat, G. Loia, J. A. Marchand

Updated: 08/22/24
"""
########################################################################
########################################################################

########################################################################
#imports
import pandas as pd
import sys
########################################################################

def load_data(file_path):
    """
    Load the TSV file into a DataFrame.
    
    Parameters:
    - file_path: str, the path to the input TSV file
    
    Returns:
    - df: pandas DataFrame, the loaded data
    """
    return pd.read_csv(file_path, sep='\t')

def tally_unique_read_positions(df):
    """
    Tally the count of unique 'read_focus_base' positions for each 'read_id'.
    
    Parameters:
    - df: pandas DataFrame, the input data
    
    Returns:
    - tally_counts: pandas Series, counts of unique positions
    """
    return df.groupby('read_id')['read_focus_base'].nunique().value_counts()

def detect_focus_range(tally_counts):
    """
    Detect the focus range based on the most frequent unique count.
    
    Parameters:
    - tally_counts: pandas Series, counts of unique positions
    
    Returns:
    - focus_range: int, the most frequent unique count
    """
    focus_range = tally_counts.idxmax()  # Most frequent count of unique positions
    #print(f'Xemora [STATUS] - Autodetected focus range: ' + str(focus_range))
    return focus_range

def filter_valid_reads(df, focus_range):
    """
    Filter out read_ids that do not have exactly 'focus_range' unique 'read_focus_base' values.
    
    Parameters:
    - df: pandas DataFrame, the input data
    - focus_range: int, the target number of unique positions
    
    Returns:
    - df_filtered: pandas DataFrame, the filtered data
    """
    read_focus_base_counts = df.groupby('read_id')['read_focus_base'].nunique()
    valid_read_ids = read_focus_base_counts[read_focus_base_counts == focus_range].index
    return df[df['read_id'].isin(valid_read_ids)]

def split_class_probs(df):
    """
    Split the 'class_probs' column into two new columns: 'prob_DNA' and 'prob_XNA'.
    
    Parameters:
    - df: pandas DataFrame, the input data
    
    Returns:
    - df: pandas DataFrame, the data with new columns added
    """
    df = df.copy()
    df[['prob_DNA', 'prob_XNA']] = df['class_probs'].str.split(',', expand=True).astype(float)
    return df

def normalize_read_focus_base(df):
    """
    Normalize 'read_focus_base' for each 'read_id', then filter out read_ids that
    do not fall within the most common range of bases analyzed.
    
    Parameters:
    - df: pandas DataFrame, the input data
    
    Returns:
    - df: pandas DataFrame, the data with normalized 'read_focus_base' and filtered read_ids
    """
    df = df.copy()

    # Step 1: Normalize the 'read_focus_base' for each 'read_id'
    df['normalized_read_focus_base'] = df.groupby('read_id')['read_focus_base'].transform(lambda x: x - x.min() + 1)

    # Step 2: Calculate the min and max of 'read_focus_base' for each 'read_id'
    ranges = df.groupby('read_id')['read_focus_base'].agg(['min', 'max'])
    
    # Step 3: Count the frequency of each (min, max) range
    range_counts = ranges.value_counts().reset_index(name='count')

    # Step 4: Find the most common range (min, max)
    most_common_range = range_counts.loc[range_counts['count'].idxmax(), ['min', 'max']]

    # Step 5: Filter read_ids that fall within the most common (min, max) range
    valid_read_ids = ranges[(ranges['min'] == most_common_range['min']) & (ranges['max'] == most_common_range['max'])].index

    # Step 6: Filter the dataframe to only include those read_ids
    df_filtered = df[df['read_id'].isin(valid_read_ids)]

    # Print the number of unique read_ids in df and df_filtered
    print(f"Number of unique read_ids in original df: {df['read_id'].nunique()}")
    print(f"Number of unique read_ids in filtered df: {df_filtered['read_id'].nunique()}")

    return df_filtered

def save_data(df, output_file):
    """
    Save the DataFrame to a CSV file.
    
    Parameters:
    - df: pandas DataFrame, the data to be saved
    - output_file: str, the path to save the CSV file
    
    Returns:
    - None
    """
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Xemora [STATUS] - Position normalization complete - per-read results saved to {output_file}")

def main():
    """
    Main function to orchestrate the data processing steps.
    """
    #Load input file name
    input_file = sys.argv[1] 
    #Set output file name
    output_file = sys.argv[2]

    # Load the input data
    df = load_data(input_file)

    # Tally unique read positions
    tally_counts = tally_unique_read_positions(df)
    #print(tally_counts)

    # Detect the focus range
    focus_range = detect_focus_range(tally_counts)

    # Filter the data based on focus range
    df_filtered = filter_valid_reads(df, focus_range)

    # Split 'class_probs' into 'prob_DNA' and 'prob_XNA'
    df_filtered = split_class_probs(df_filtered)

    # Normalize 'read_focus_base'
    df_filtered = normalize_read_focus_base(df_filtered)

    # Save the processed data to a CSV file
    save_data(df_filtered, output_file)

    # Print the filtered and normalized DataFrame
    #print(df_filtered)

    # Re-tally the counts after filtering
    final_tally_counts = tally_unique_read_positions(df_filtered)
    #print(final_tally_counts)

if __name__ == "__main__":
    main()

