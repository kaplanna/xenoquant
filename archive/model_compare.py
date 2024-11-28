import pandas as pd
import os

def compare_model_outputs(GC_Model, AT_Model):
    """
    Compare the outputs of two models and save the comparison results.
    For each alignment_results file:
        - Loads the data
        - Merges based on read_id
        - Applies AND gate logic for final prediction
        - Saves the comparison results
        - Counts and prints the number and fractions of reads categorized as XNA, DNA, AT, GC
    """
    # Ensure the output directories exist or create them
    output_dir_GC = os.path.join(GC_Model, 'model_compare_results')
    output_dir_AT = os.path.join(AT_Model, 'model_compare_results')
    os.makedirs(output_dir_GC, exist_ok=True)
    os.makedirs(output_dir_AT, exist_ok=True)
    
    overall_results = []

    for file_name in os.listdir(GC_Model):
        if file_name.startswith('alignment_results'):
            # Load alignment files from two different directories
            file_path_GC = os.path.join(GC_Model, file_name)
            file_path_AT = os.path.join(AT_Model, file_name)
            
            # Load the data into dataframes
            df_GC = pd.read_csv(file_path_GC)
            df_AT = pd.read_csv(file_path_AT)
            
            # Merge the two dataframes based on read_id
            merged_df = pd.merge(df_GC, df_AT, on='read_id', suffixes=('_GC', '_AT'))
            
            # Apply AND gate logic for final prediction
            merged_df['final_prediction'] = merged_df['class_pred_GC'] & merged_df['class_pred_AT']
            
            # Select relevant columns to include in the output
            output_df = merged_df[['read_id', 'class_pred_GC', 'class_pred_AT', 'class_0_probs_GC', 'class_1_probs_GC', 'class_0_probs_AT', 'class_1_probs_AT', 'final_prediction']]
            
            # Save the results in the new directories
            output_df.to_csv(os.path.join(output_dir_GC, f'{file_name}_comparison.csv'), index=False)
            output_df.to_csv(os.path.join(output_dir_AT, f'{file_name}_comparison.csv'), index=False)
            
            # Count the reads and calculate fractions for different categories
            total_calls = len(output_df)
            xna_count = output_df['final_prediction'].sum()
            dna_count = total_calls - xna_count
            at_count = len(output_df[output_df['class_pred_AT'] == 0])
            gc_count = len(output_df[output_df['class_pred_GC'] == 0])
            
            # Calculate fractions
            xna_frac = xna_count / total_calls
            dna_frac = dna_count / total_calls
            at_frac = at_count / total_calls
            gc_frac = gc_count / total_calls
            # Save individual file results and append data to the overall results
            overall_results.append({
                'file_name': file_name,
                'XNA_count': xna_count, 'XNA_fraction': xna_frac,
                'DNA_count': dna_count, 'DNA_fraction': dna_frac,
                'AT_count': at_count, 'AT_fraction': at_frac,
                'GC_count': gc_count, 'GC_fraction': gc_frac,
            })
            
            
            # Print the counts and fractions for each category
            print(f"For file {file_name}:")
            print(f"Number of reads called as PZ: {xna_count} (Fraction: {xna_frac:.2f})")
            print(f"Number of reads called as GC: {gc_count} (Fraction: {gc_frac:.2f})")
            print(f"Number of reads called as AT: {at_count} (Fraction: {at_frac:.2f})\n")
            
            
            
    # Save the overall results as a CSV file in each output directory
    overall_df = pd.DataFrame(overall_results)
    overall_df.to_csv(os.path.join(output_dir_GC, 'overall_comparison_results.csv'), index=False)
    overall_df.to_csv(os.path.join(output_dir_AT, 'overall_comparison_results.csv'), index=False)

# Specify the results directories
GC_Model = '/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/XPCR_GAP_Basecall/ZC_Basecall/Results'
AT_Model = '/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/XPCR_GAP_Basecall/ZT_Basecall/Results'


compare_model_outputs(GC_Model, AT_Model)

