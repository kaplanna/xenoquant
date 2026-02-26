import pysam
import pandas as pd

def export_fasta_from_bam(bam_path, read_id_csv, output_fasta):
    # Load the read ID pairs
    df = pd.read_csv(read_id_csv)  # Change to ',' if comma-separated

    # Flatten both columns into one list of unique read IDs
    read_ids = pd.unique(df[['read_id_1', 'read_id_2']].values.ravel())

    # Open the BAM file
    bam = pysam.AlignmentFile(bam_path, "rb")

    # Open the FASTA output file
    with open(output_fasta, "w") as fasta_out:
        for read in bam:
            if read.query_name in read_ids and not read.is_unmapped:
                sequence = read.query_sequence
                fasta_out.write(f">{read.query_name}\n{sequence}\n")

    print(f"FASTA file saved to: {output_fasta}")

# Example usage
export_fasta_from_bam(
    bam_path="/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/XPCR/250327_B45_Duplex/duplex_output/output.bam",
    read_id_csv="/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/XPCR/250327_B45_Duplex/duplex_output/duplex_valid_only.csv",
    output_fasta="/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/XPCR/250327_B45_Duplex/duplex_output/duplex_sequences.fasta"
)

