import pandas as pd
import pysam

def get_read_lengths(bam_path, read_id_pairs_csv, output_csv):
    # Load the read ID pairs
    df = pd.read_csv(read_id_pairs_csv)  # Change to ',' if your CSV is comma-delimited

    # Open the BAM file
    bam = pysam.AlignmentFile(bam_path, "rb")

    # Create dictionaries to store read lengths
    read_lengths = {}

    # Iterate over all reads in the BAM file and store lengths
    for read in bam:
        if not read.is_unmapped:
            read_lengths[read.query_name] = read.query_length

    # Map read lengths to the input DataFrame
    df['read_1_length'] = df['read_id_1'].map(read_lengths)
    df['read_2_length'] = df['read_id_2'].map(read_lengths)

    # Save the output
    df.to_csv(output_csv, index=False)
    print(f"Saved read lengths to {output_csv}")

# Example usage
get_read_lengths(
    bam_path='/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/XPCR/250327_B45_Duplex/duplex_output/output.bam',
    read_id_pairs_csv='/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/XPCR/250327_B45_Duplex/duplex_output/duplex_valid_only.csv',
    output_csv='/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/XPCR/250327_B45_Duplex/duplex_output/read_lengths_output.csv'
)


