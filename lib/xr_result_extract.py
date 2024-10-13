import os
import pandas as pd

def combine_results(parent_dir, output_file):
    combined_df = pd.DataFrame()  # Empty DataFrame to hold all results
    
    # Iterate through each folder in the parent directory
    for folder_name in os.listdir(parent_dir):
        folder_path = os.path.join(parent_dir, folder_name)
        
        # Check if it's a directory
        if os.path.isdir(folder_path):
            results_file = os.path.join(folder_path, 'demux', 'overall_demux_results.csv')
            
            # Check if the file exists
            if os.path.exists(results_file):
                # Read the CSV file
                df = pd.read_csv(results_file)
                
                # Add a column with the folder name
                df['folder_name'] = folder_name
                
                # Append to the combined DataFrame
                combined_df = pd.concat([combined_df, df], ignore_index=True)
    
    # Write the combined DataFrame to a new CSV file
    output_file = os.path.join(parent_dir, output_csv)
    combined_df.to_csv(output_file, index=False)
    print(f"Combined CSV saved to {output_file}")

# Example usage
parent_directory = '/home/xenolab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/241003_B5R_Basecall'  
output_csv = 'full_results.csv'
combine_results(parent_directory, output_csv)

