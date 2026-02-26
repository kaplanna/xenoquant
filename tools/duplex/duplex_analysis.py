import pysam
from collections import defaultdict

def extract_related_read_ids(bam_file, output_file):
    """
    Extracts read IDs that are related by duplex (based on the `dx` tag).
    
    Parameters:
        bam_file (str): Path to the input BAM file.
        output_file (str): Path to the output text file containing grouped read IDs.
    """
    # Open the BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bam:
    # Process BAM file

        related_reads = defaultdict(list)

        # Iterate over each read in the BAM file
        for read in bam:
            if read.has_tag("dx"):
                dx_value = read.get_tag("dx")
                read_id = read.query_name
                related_reads[dx_value].append(read_id)

    # Write the results to the output file
    with open(output_file, "w") as f:
        for dx_value, read_ids in related_reads.items():
            f.write(f"dx:{dx_value}\n")
            f.write("\n".join(read_ids) + "\n\n")
    
    print(f"Related read IDs written to {output_file}")

# Example usage
if __name__ == "__main__":
    bam_file = "/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/XPCR/250327_B45_Duplex/duplex_output/output.bam"  # Replace with your BAM file path
    output_file = "/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/XPCR/250327_B45_Duplex/duplex_output/duplex_read_IDs.txt"  # Replace with your desired output file path
    
    extract_related_read_ids(bam_file, output_file)

