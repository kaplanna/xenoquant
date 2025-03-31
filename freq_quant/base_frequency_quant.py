import pysam
import pandas as pd
from collections import Counter
import os

# Define file paths directly in the script
BAM_FILE = "/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/XPCR/250315_B43-B44_Basecall/GBC_BA-Basecall/preprocess/bam/aligned.BAM"
BED_FILE = "/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/XPCR/250315_B43-B44_Basecall/GBC_BA-Basecall/references/B.bed"
OUTPUT_CSV = "/home/marchandlab/github/kaplanna/xemora/freq_quant/B43_B44_output_analysis.csv"
ALL_READ_IDS_CSV = "/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/XPCR/250315_B43-B44_Basecall/GBC_BA-Basecall/demux/all_read_ids.csv"

def ensure_bam_index(bam_file):
    """Ensure that the BAM file has an index; create one if missing."""
    index_file = bam_file + ".bai"
    if not os.path.exists(index_file):
        print(f"Index file not found for {bam_file}. Generating index...")
        pysam.index(bam_file)
    else:
        print(f"Index file found: {index_file}")

def extract_positions_from_bed(bed_file):
    """Extract reference positions from a BED file."""
    positions = []
    with open(bed_file, 'r') as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) >= 3:
                ref_name = cols[0]
                position = int(cols[1]) + 1  # BED is 0-based, convert to 1-based
                positions.append((ref_name, position))
    return positions

def analyze_position(bam_file, reference_name, position):
    """
    Analyze base and indel frequencies at a given position in an aligned BAM file.
    
    Parameters:
        bam_file (str): Path to the BAM file.
        reference_name (str): Reference sequence name (e.g., chromosome or contig).
        position (int): 1-based reference position to analyze.
    
    Returns:
        list: A list of dictionaries containing read_id, base, and barcode.
    """
    bam = pysam.AlignmentFile(bam_file, "rb")
    results = []
    
    for read in bam.fetch(reference_name, position - 1, position):
        read_id = read.query_name
        for query_pos, ref_pos in read.get_aligned_pairs(matches_only=False):
            if ref_pos == position - 1:
                base = read.query_sequence[query_pos] if query_pos is not None else 'DEL'
                results.append({
                    'read_id': read_id,
                    'base': base
                })
        
        # Check for insertions
        for cigartuples in read.cigartuples:
            if cigartuples[0] == 1:  # 1 = Insertion
                results.append({
                    'read_id': read_id,
                    'base': 'INS'
                })
    
    bam.close()
    return results

def main():
    ensure_bam_index(BAM_FILE)
    positions = extract_positions_from_bed(BED_FILE)
    
    results = []
    for ref_name, position in positions:
        results.extend(analyze_position(BAM_FILE, ref_name, position))
    
    # Convert to DataFrame
    df = pd.DataFrame(results)
    
    # Load the barcode mapping file
    barcode_df = pd.read_csv(ALL_READ_IDS_CSV, sep=',')
    if 'barcode_pair,read_id' in barcode_df.columns:
        barcode_df[['barcode_pair', 'read_id']] = barcode_df['barcode_pair,read_id'].str.split(',', expand=True)
        barcode_df.drop(columns=['barcode_pair,read_id'], inplace=True)
    
    # Merge on read_id
    df = df.merge(barcode_df, on='read_id', how='left')
    
    # Group by barcode and count bases
    summary_df = df.groupby('barcode_pair')['base'].value_counts().unstack(fill_value=0)
    
    # Reorder columns
    desired_order = ['INS', 'DEL', 'A', 'T', 'G', 'C']
    summary_df = summary_df.reindex(columns=desired_order, fill_value=0)
    
    # Save to CSV
    summary_df.to_csv(OUTPUT_CSV)
    print(f"Summary results saved to {OUTPUT_CSV}")

if __name__ == "__main__":
    main()

