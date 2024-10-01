import os
import subprocess

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in reversed(seq)])

def run_cutadapt(input_fastq, adapter, output_file, error_rate, min_overlap, min_len, max_len):
    cmd = [
        'cutadapt',
        '-g', adapter,  # Forward barcode
        '-m', str(min_len),  # Minimum length
        '-M', str(max_len),  # Maximum length
        '--no-trim',  # Do not trim off the barcode sequence
        '--discard-untrimmed',  # Discard untrimmed reads
        '-e', str(error_rate),  # Error rate
        '-O', str(min_overlap),  # Minimum overlap
        '-o', output_file,  # Output file
        input_fastq  # Input FASTQ or bam file
    ]
    subprocess.run(cmd, check=True)

def combine_and_deduplicate_fastq(file1, file2, output_file):
    unique_reads = set()
    with open(file1, 'r') as infile1, open(file2, 'r') as infile2, open(output_file, 'w') as outfile:
        for infile in [infile1, infile2]:
            for title, seq, plus, qual in zip(*[infile]*4):
                if seq not in unique_reads:
                    unique_reads.add(seq)
                    outfile.writelines([title, seq, plus, qual])

def parse_barcodes(barcode_file):
    forward_barcodes = {}
    reverse_barcodes = {}
    with open(barcode_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                name = line.strip()[1:]
                sequence = next(f).strip()
                if 'FWD' in name:
                    forward_barcodes[name] = sequence
                elif 'REV' in name:
                    reverse_barcodes[name] = sequence
    return forward_barcodes, reverse_barcodes

def generate_barcode_pairs(forward_barcodes, reverse_barcodes):
    pairs = {}
    for fwd_name, fwd_seq in forward_barcodes.items():
        for rev_name, rev_seq in reverse_barcodes.items():
            pair_name = f"{fwd_name}_{rev_name}"
            pairs[pair_name] = (fwd_seq, rev_seq)
    return pairs

# Define the paths to the input files and output directory
input_fastq = "/home/xenolab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/240930_NTC_Phusion_Training_Testing/240930_NTC_Phusion_Withheld_Test_Set/bc.bam"
output_dir = "/home/xenolab/DataAnalysis/Kaplan/basecall/10.4.1/BSn/240930_NTC_Phusion_Training_Testing/240930_NTC_Phusion_Withheld_Test_Set/demux"
barcode_file = "/home/xenolab/github/kaplanna/xemora/training_demux/NTC_barcodes_CA.fasta"

os.makedirs(output_dir, exist_ok=True)

error_rate = 0.15
min_overlap = 18
min_len = 50
max_len = 250

# Parse the barcodes from the file
forward_barcodes, reverse_barcodes = parse_barcodes(barcode_file)

# Generate all forward-reverse barcode pairs
barcode_pairs = generate_barcode_pairs(forward_barcodes, reverse_barcodes)

for barcode_pair, (barcode1, barcode2) in barcode_pairs.items():
    barcode2_rc = reverse_complement(barcode2)
    barcode1_rc = reverse_complement(barcode1)
    combined_adapter1 = f'{barcode1}...{barcode2_rc}'
    combined_adapter2 = f'{barcode2}...{barcode1_rc}'

    output_file1 = os.path.join(output_dir, f'{barcode_pair}_1.fastq')
    output_file2 = os.path.join(output_dir, f'{barcode_pair}_2.fastq')

    print(f'Running cutadapt for barcode pair: {barcode_pair}, adapter 1')
    run_cutadapt(input_fastq, combined_adapter1, output_file1, error_rate, min_overlap, min_len, max_len)

    print(f'Running cutadapt for barcode pair: {barcode_pair}, adapter 2')
    run_cutadapt(input_fastq, combined_adapter2, output_file2, error_rate, min_overlap, min_len, max_len)

    # Combine and deduplicate outputs
    dedup_output = os.path.join(output_dir, f'{barcode_pair}.fastq')
    print(f'Combining and deduplicating output for barcode pair: {barcode_pair}')
    combine_and_deduplicate_fastq(output_file1, output_file2, dedup_output)

    # Delete the non-deduplicated files
    os.remove(output_file1)
    os.remove(output_file2)

print('Cutadapt processing, deduplication, and cleanup completed.')


