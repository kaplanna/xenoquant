
import os
from Bio import SeqIO

def extract_read_ids_per_fastq(working_dir):
    """
    Loops through FASTQ files in a directory and writes a line-delimited text file for each FASTQ file with all read IDs.
    Also prints and logs the number of reads in each FASTQ file.

    Parameters:
    working_dir (str): Path to the directory containing FASTQ files.

    Returns:
    None
    """
    
    # Create a demux summary file to store the number of reads per FASTQ file
    demux_summary_file = os.path.join(working_dir, "demux.txt")
    
    with open(demux_summary_file, 'w') as summary_file:
        summary_file.write("FASTQ File\tNumber of Reads\n")
        
        # Loop through each FASTQ file in the directory
        for fastq_file in os.listdir(working_dir):
            if fastq_file.endswith('.fastq') or fastq_file.endswith('.fq'):
                fastq_path = os.path.join(working_dir, fastq_file)
                output_txt_file = os.path.join(working_dir, f"{os.path.splitext(fastq_file)[0]}_read_ids.txt")
                print(f"Processing {fastq_path}...")

                read_count = 0
                # Open the FASTQ file and extract read IDs
                try:
                    with open(fastq_path, 'r') as fq, open(output_txt_file, 'w') as out_file:
                        for record in SeqIO.parse(fq, "fastq"):
                            read_id = record.id
                            out_file.write(f"{read_id}\n")
                            read_count += 1

                    # Print and write the read count to the summary file
                    print(f"{fastq_file} contains {read_count} reads.")
                    summary_file.write(f"{fastq_file}\t{read_count}\n")
                    print(f"Read IDs for {fastq_file} written to {output_txt_file}")

                except Exception as e:
                    print(f"Error processing {fastq_path}: {e}")

working_dir = "/home/xenolab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/240930_NTC_Phusion_Training_Testing/240930_NTC_Phusion_750_Training_Set/demux"
extract_read_ids_per_fastq(working_dir)
