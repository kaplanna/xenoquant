import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.stats import ttest_ind
from adjustText import adjust_text

def normalize_frac(df, nc_sample='C-AT', pc_sample='C-BS'):
    """
    Normalizes XNA fractions for B and S between the same positive and negative controls.

    Parameters:
    df (pd.DataFrame): The input DataFrame.
    nc_sample (str): The sample identifier for the negative control for both B and S (default 'C-AT').
    pc_sample (str): The sample identifier for the positive control for both B and S (default 'C-BS').

    Returns:
    pd.DataFrame: The modified DataFrame with additional columns for normalized XNA fractions B and S.
    """
    # Identify the XNA fractions for the controls for B and S
    nc_value_b = df.loc[df['Sample'] == nc_sample, 'XNA_fraction_B'].values[0]
    pc_value_b = df.loc[df['Sample'] == pc_sample, 'XNA_fraction_B'].values[0]

    nc_value_s = df.loc[df['Sample'] == nc_sample, 'XNA_fraction_S'].values[0]
    pc_value_s = df.loc[df['Sample'] == pc_sample, 'XNA_fraction_S'].values[0]

    # Normalize the XNA fractions for B and S between the PC and NC
    df['Normalized_XNA_fraction_B'] = (df['XNA_fraction_B'] - nc_value_b) / (pc_value_b - nc_value_b)
    df['Normalized_XNA_fraction_S'] = (df['XNA_fraction_S'] - nc_value_s) / (pc_value_s - nc_value_s)
    
    return df

def load_and_merge_data(frac_file, map_file):
    """
    Loads fraction and mapping data, splits Sample ID, and merges data frames on barcodes.

    Parameters:
    frac_file (str): The path to the fraction CSV file.
    map_file (str): The path to the PCR mapping CSV file.

    Returns:
    pd.DataFrame: The merged DataFrame.
    """
    # Load fraction file as DataFrame
    frac_df = pd.read_csv(frac_file)

    # Load PCR mapping file as DataFrame
    map_df = pd.read_csv(map_file)
    map_df[['Experiment', 'Sample', 'Replicate']] = map_df['Sample ID'].str.split('.', expand=True)

    # Merge the DataFrames on the FWD and REV barcodes
    merged_df = pd.merge(map_df, frac_df, left_on=['BC-FWD', 'BC-REV'], right_on=['FWD_barcode', 'REV_barcode'])

    return merged_df

def get_consistent_columns(df):
    """
    Identifies consistent columns across replicates within each group of Experiment and Sample.

    Parameters:
    df (pd.DataFrame): The DataFrame to check for consistency.

    Returns:
    list: A list of consistent column names.
    """
    # Group by Experiment and Sample
    grouped = df.groupby(['Experiment', 'Sample'])

    # Identify columns that have only one unique value within each group
    consistent_columns = []
    for col in df.columns:
        if grouped[col].nunique().max() == 1:
            consistent_columns.append(col)
    
    # Ensure essential columns are included
    essential_columns = ['Experiment', 'Sample']
    consistent_columns = list(set(consistent_columns) | set(essential_columns))

    print(f"Consistent columns: {consistent_columns}")
    return consistent_columns

def avg_stdev_calc(df, consistent_columns):
    """
    Calculates the average and standard deviation for each sample across replicates,
    for both normalized and non-normalized XNA fractions, keeping consistent columns.

    Parameters:
    df (pd.DataFrame): The DataFrame containing normalized and non-normalized XNA fractions.
    consistent_columns (list): The list of consistent column names.

    Returns:
    pd.DataFrame: A summary DataFrame with averages and standard deviations, sorted by Experiment and Sample.
    """
    # Group by consistent columns and calculate mean and std for both normalized and non-normalized XNA fractions
    summary_df = df.groupby(consistent_columns).agg(
        avg_xna_b=pd.NamedAgg(column='XNA_fraction_B', aggfunc='mean'),
        std_xna_b=pd.NamedAgg(column='XNA_fraction_B', aggfunc='std'),
        avg_xna_s=pd.NamedAgg(column='XNA_fraction_S', aggfunc='mean'),
        std_xna_s=pd.NamedAgg(column='XNA_fraction_S', aggfunc='std'),
        avg_normalized_xna_b=pd.NamedAgg(column='Normalized_XNA_fraction_B', aggfunc='mean'),
        std_normalized_xna_b=pd.NamedAgg(column='Normalized_XNA_fraction_B', aggfunc='std'),
        avg_normalized_xna_s=pd.NamedAgg(column='Normalized_XNA_fraction_S', aggfunc='mean'),
        std_normalized_xna_s=pd.NamedAgg(column='Normalized_XNA_fraction_S', aggfunc='std')
    ).reset_index()

    # Sort by Experiment and Sample
    summary_df = summary_df.sort_values(by=['Experiment', 'Sample'])

    return summary_df

def format_data_for_volcano(df, output_file='volcano_data.csv'):
    experiments = df['Experiment'].unique()
    condition_list = []
    model_list = []
    experiment_list = []
    rep1_list = []
    rep2_list = []
    rep3_list = []
    ref1_list = []
    ref2_list = []
    ref3_list = []

    # Filter out control and reference samples from the condition list
    df_conditions = df[df['Type'] == 'Sample']
    df_references = df[df['Type'] == 'Reference']

    for experiment in experiments:
        exp_df = df_conditions[df_conditions['Experiment'] == experiment]
        ref_df = df_references[df_references['Experiment'] == experiment]

        conditions = exp_df['Sample'].unique()

        for condition in conditions:
            cond_df = exp_df[exp_df['Sample'] == condition]

            for model, fraction_column in [('BA', 'Normalized_XNA_fraction_B'), ('ST', 'Normalized_XNA_fraction_S')]:
                # Extract replicates for the condition
                replicates = cond_df.sort_values('Replicate')[[fraction_column]].values
                references = ref_df[[fraction_column]].values

                # Ensure we have exactly 3 replicates
                replicates = np.vstack([replicates, np.full((3 - len(replicates), 1), np.nan)])[:3]
                references = np.vstack([references, np.full((3 - len(references), 1), np.nan)])[:3]

                # Append the data for this condition and model
                condition_list.append(f"{condition}")
                model_list.append(model)
                experiment_list.append(experiment)
                rep1_list.append(replicates[0][0])
                rep2_list.append(replicates[1][0])
                rep3_list.append(replicates[2][0])
                ref1_list.append(references[0][0])
                ref2_list.append(references[1][0])
                ref3_list.append(references[2][0])

    # Create the DataFrame
    formatted_df = pd.DataFrame({
        'Condition': condition_list,
        'Model': model_list,
        'Experiment': experiment_list,
        'Rep1': rep1_list,
        'Rep2': rep2_list,
        'Rep3': rep3_list,
        'Ref1': ref1_list,
        'Ref2': ref2_list,
        'Ref3': ref3_list,
    })

    # Save to CSV
    formatted_df.to_csv(output_file, index=False)
    print(f"Formatted data saved to {output_file}")


def plot_bar_plots(df, output_dir='./plots'):
    """
    Generates bar plots for each condition with both normalized B and S XNA fractions on the same plot, with standard deviation error bars.

    Parameters:
    df (pd.DataFrame): The summary DataFrame containing averages and standard deviations.
    output_dir (str): The directory to save the plots (default './plots').
    """
    # Create output directory if it doesn't exist
    import os
    os.makedirs(output_dir, exist_ok=True)

    # Set up the figure and axes
    fig, ax = plt.subplots(figsize=(12, 8))

    # Set the bar width
    bar_width = 0.35

    # Calculate positions for the bars
    indices = np.arange(len(df))

    # Plot bars for B
    bars_b = ax.bar(indices, df['avg_normalized_xna_b'], bar_width, yerr=df['std_normalized_xna_b'], label='Normalized XNA B', color='skyblue', capsize=4)

    # Plot bars for S, with an offset
    bars_s = ax.bar(indices + bar_width, df['avg_normalized_xna_s'], bar_width, yerr=df['std_normalized_xna_s'], label='Normalized XNA S', color='lightgreen', capsize=4)

    # Add labels and title
    ax.set_xlabel('Sample')
    ax.set_ylabel('Normalized XNA Fraction')
    ax.set_title('Normalized XNA Fractions B and S by Sample')
    ax.set_xticks(indices + bar_width / 2)
    ax.set_xticklabels(df['Sample'], rotation=45, ha='right')

    # Add a legend
    ax.legend()

    # Adjust layout to make room for the x-axis labels
    plt.tight_layout()

    # Save the plot
    plt.savefig(os.path.join(output_dir, 'normalized_xna_b_s.png'))

    # Show the plot
    plt.show()

def plot_sample_vs_reference(df, output_dir='./plots'):
    """
    Compares each sample to its experimental reference condition and generates bar plots.

    Parameters:
    df (pd.DataFrame): The DataFrame containing normalized XNA fractions and reference information.
    output_dir (str): The directory to save the plots (default './plots').
    """
    # Create output directory if it doesn't exist
    import os
    os.makedirs(output_dir, exist_ok=True)

    # Filter reference conditions
    reference_df = df[df['Type'] == 'Reference'][['Experiment', 'avg_normalized_xna_b', 'avg_normalized_xna_s']]
    reference_df = reference_df.rename(columns={'avg_normalized_xna_b': 'ref_normalized_xna_b', 'avg_normalized_xna_s': 'ref_normalized_xna_s'})

    # Merge with original data to get corresponding references
    comparison_df = pd.merge(df[df['Type'] == 'Sample'], reference_df, on='Experiment')

    # Calculate differences
    comparison_df['diff_normalized_xna_b'] = comparison_df['avg_normalized_xna_b'] - comparison_df['ref_normalized_xna_b']
    comparison_df['diff_normalized_xna_s'] = comparison_df['avg_normalized_xna_s'] - comparison_df['ref_normalized_xna_s']

    # Plotting
    fig, ax = plt.subplots(figsize=(14, 8))

    # Set the bar width
    bar_width = 0.35

    # Calculate positions for the bars
    indices = np.arange(len(comparison_df))

    # Plot bars for B difference
    bars_b = ax.bar(indices, comparison_df['diff_normalized_xna_b'], bar_width, label='Difference B', color='blue', capsize=4)

    # Plot bars for S difference, with an offset
    bars_s = ax.bar(indices + bar_width, comparison_df['diff_normalized_xna_s'], bar_width, label='Difference S', color='green', capsize=4)

    # Add labels and title
    ax.set_xlabel('Sample')
    ax.set_ylabel('Difference from Reference')
    ax.set_title('Difference in Normalized XNA Fractions B and S Compared to Reference')
    ax.set_xticks(indices + bar_width / 2)
    ax.set_xticklabels(comparison_df['Sample'], rotation=45, ha='right')

    # Add a legend
    ax.legend()

    # Adjust layout to make room for the x-axis labels
    plt.tight_layout()

    # Save the plot
    plt.savefig(os.path.join(output_dir, 'sample_vs_reference.png'))

    # Show the plot
    plt.show()




def main():
    # File paths
    frac_file = '/Users/nickkaplan/DataAnalysis/basecall/240830_B16_B17_Basecall/demux/output_xna_fractions.csv'
    map_file = '/Users/nickkaplan/DataAnalysis/basecall/240830_B16_B17_Basecall/B16_B17_map.csv'
    output_file = '/Users/nickkaplan/DataAnalysis/basecall/240830_B16_B17_Basecall/normalized_xna_fraction_output.csv'
    summary_file = '/Users/nickkaplan/DataAnalysis/basecall/240830_B16_B17_Basecall/summary_xna_fractions.csv'
    volcano_file = '/Users/nickkaplan/DataAnalysis/basecall/240830_B16_B17_Basecall/volcano_data.csv'


    # Load and merge the data
    merged_df = load_and_merge_data(frac_file, map_file)

    # Normalize the XNA fractions
    normalized_df = normalize_frac(merged_df)

    # Identify consistent columns across replicates
    consistent_columns = get_consistent_columns(normalized_df)

    # Save the final cleaned and normalized DataFrame
    normalized_df.to_csv(output_file, index=False)
    print(f"Normalized and cleaned file saved as {output_file}")
    
    # Format and save the volcano data
    format_data_for_volcano(normalized_df, output_file=volcano_file)

    # Create summary with average and standard deviation, using only consistent columns
    summary_df = avg_stdev_calc(normalized_df, consistent_columns)
    summary_df.to_csv(summary_file, index=False)
    print(f"Summary file saved as {summary_file}")

    # Generate plots
    plot_bar_plots(summary_df)
    plot_sample_vs_reference(summary_df)

    print("Plots saved successfully.")

if __name__ == "__main__":
    main()
