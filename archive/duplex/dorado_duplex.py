import os
import subprocess

def run_dorado_duplex(dorado_path, pod5_file, output_bam, reference_fasta):
    """
    Runs the dorado duplex command on a specified .pod5 file using a specific dorado executable.
    
    Parameters:
        dorado_path (str): Path to the dorado executable.
        pod5_file (str): Path to the input .pod5 file.
        output_bam (str): Path to the output BAM file.
    """
    # Check if the dorado executable exists
    if not os.path.isfile(dorado_path):
        print(f"Error: The specified dorado executable does not exist: {dorado_path}")
        return
    
    # Check if the input pod5 file exists
    if not os.path.isfile(pod5_file):
        print(f"Error: The specified .pod5 file does not exist: {pod5_file}")
        return
    
    # Build the dorado duplex command
    command = f"{dorado_path} duplex sup --reference {reference_fasta} {pod5_file} > {output_bam}"
    print(f"Running command: {command}")
    try:
        # Execute the command
        subprocess.run(command, shell=True, check=True)
        print(f"Duplex basecalling completed. Output saved to: {output_bam}")
    except subprocess.CalledProcessError as e:
        print(f"Error running dorado duplex: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# Example usage
if __name__ == "__main__":
    # Specify the path to the dorado executable, input pod5 file, and output BAM file
    dorado_path = "/home/marchandlab/dorado-0.8.0-linux-x64/bin/dorado"  # Replace with the path to your dorado executable
    pod5_file = "/home/marchandlab/DataAnalysis/Kaplan/raw/xPCR/250326_B45_OEP_Duplex/20250326_1806_MN37138_AYF750_2737b5b8/merged_pod5/pod5.pod5"  # Replace with your .pod5 file path
    output_bam = "/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/XPCR/250327_B45_Duplex/duplex_output/output.bam"  # Replace with desired output BAM file path
    reference_fasta = "/home/marchandlab/DataAnalysis/Kaplan/raw/xPCR/250326_B45_OEP_Duplex/ref/ref_B45.fasta"

    # Run the dorado duplex command
    run_dorado_duplex(dorado_path, pod5_file, output_bam, reference_fasta)

