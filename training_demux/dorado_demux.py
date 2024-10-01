import subprocess
import os
import pysam

def dorado_demux(dorado_path, working_dir, barcode_arr_file, barcode_file, bam_path):
    """
    Function to demultiplex a basecall file using Dorado demux with custom barcode files.

    Parameters:
    dorado_path (str): Path to the Dorado executable.
    working_dir (str): Path to the working directory for output.
    barcode_arr_file (str): Path to the TOML file with barcode arrangement configuration.
    barcode_file (str): Path to the FASTA file containing the barcode sequences.
    bam_path (str): Path to the basecall file (FAST5 or BAM format).

    Returns:
    None
    """
    
    # Ensure working directory exists
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)

    # Construct the dorado demux command
    barcode_cmd = f"{dorado_path} demux --output-dir {working_dir} --barcode-arrangement {barcode_arr_file} --barcode-sequences {barcode_file} {bam_path}"

    try:
        # Run the command
        subprocess.run(barcode_cmd, shell=True, check=True)
        print(f"Demultiplexing completed successfully. Output saved in {working_dir}")
    except subprocess.CalledProcessError as e:
        print(f"Error during demultiplexing: {e}")


def extract_read_ids_per_bam(working_dir):
    """
    Loops through BAM files in a directory and writes a line-delimited text file for each BAM file with all read IDs.
    Also prints and logs the number of reads in each BAM file.

    Parameters:
    working_dir (str): Path to the directory containing BAM files.

    Returns:
    None
    """
    
    # Create a demux summary file to store the number of reads per BAM file
    demux_summary_file = os.path.join(working_dir, "demux.txt")
    
    with open(demux_summary_file, 'w') as summary_file:
        summary_file.write("BAM File\tNumber of Reads\n")
        
        # Loop through each BAM file in the directory
        for bam_file in os.listdir(working_dir):
            if bam_file.endswith('.bam'):
                bam_path = os.path.join(working_dir, bam_file)
                output_txt_file = os.path.join(working_dir, f"{os.path.splitext(bam_file)[0]}_read_ids.txt")
                print(f"Processing {bam_path}...")

                read_count = 0
                # Open the BAM file using pysam with check_sq=False to avoid the missing reference error
                try:
                    with pysam.AlignmentFile(bam_path, 'rb', check_sq=False) as bam, open(output_txt_file, 'w') as out_file:
                        for read in bam:
                            read_id = read.query_name
                            out_file.write(f"{read_id}\n")
                            read_count += 1

                    # Print and write the read count to the summary file
                    print(f"{bam_file} contains {read_count} reads.")
                    summary_file.write(f"{bam_file}\t{read_count}\n")
                    print(f"Read IDs for {bam_file} written to {output_txt_file}")

                except ValueError as e:
                    print(f"Error processing {bam_path}: {e}")


# Example usage
dorado_path = '~/dorado-0.7.2-linux-x64/bin/dorado'
working_dir = "/home/xenolab/DataAnalysis/Kaplan/basecall/240627_NTC_Phusion_xr_Train_Basecall/demux"
barcode_arr_file = "/home/xenolab/DataAnalysis/Kaplan/demux/barcodes.toml"
barcode_file = "/home/xenolab/DataAnalysis/Kaplan/demux/NTC_barcodes.fasta"
bam_path = "/home/xenolab/DataAnalysis/Kaplan/basecall/240627_NTC_Phusion_xr_Train_Basecall/bc.bam"

# Run Dorado demux
dorado_demux(dorado_path, working_dir, barcode_arr_file, barcode_file, bam_path)

# Extract read IDs and generate summary
extract_read_ids_per_bam(working_dir)

