import os
import numpy as np
import json
import pysam
import pod5
from remora import io
import pandas as pd

# --- Set working directory ---
working_dir = '/home/marchandlab/DataAnalysis/Kaplan/basecall/10.4.1/Manuscript_Model_Testing/BS/CBG/Can/BA-Basecall'

# --- Configurable parameters ---
MAX_READS = 100
XNA = 'B'
STRAND_FILTER = "-"  # Options: "+", "-", "both"


# --- Define fixed input paths ---
bed = os.path.join(working_dir, 'references', f'{XNA}.bed')



preprocess_dir = os.path.join(working_dir, 'preprocess')
bam_path = os.path.join(preprocess_dir, 'bam', 'aligned.BAM')
pod5_dir = os.path.join(preprocess_dir, 'pod5')


# --- Load XNA positions independently ---
xna_pos = int(pd.read_csv(bed, sep="\t", header=None).iloc[0, 1])

print(f"XNA position: {xna_pos}")


# Find the one POD5 file in the pod5 directory
pod5_files = [f for f in os.listdir(pod5_dir) if f.endswith('.pod5')]
if len(pod5_files) != 1:
    raise RuntimeError(f"Expected 1 POD5 file, found {len(pod5_files)} in {pod5_dir}")
pod5_path = os.path.join(pod5_dir, pod5_files[0])

# Create output folder and file path
output_dir = os.path.join(working_dir, 'vis_extract')
os.makedirs(output_dir, exist_ok=True)
output_json = os.path.join(output_dir, 'B_mod_output.json')


# --- Read and parse files ---
bam_file = pysam.AlignmentFile(bam_path, "rb")
pod5_dr = pod5.Reader(pod5_path)
data = []

np.set_printoptions(threshold=np.inf)

def clean_array_to_string(array):
    return ','.join(map(str, array))

# --- Extract data ---
for i, bam_read in enumerate(bam_file):
    if len(data) >= MAX_READS:
        print(f"⏹️ Reached {MAX_READS} processed reads. Stopping early.")
        break

    if bam_read.is_unmapped or bam_read.is_secondary or bam_read.is_supplementary:
        continue
    if bam_read.mapping_quality < 60:
        continue  # or a threshold like 30 for more inclusiveness
   


    read_id = bam_read.query_name
    try:
        pod5_read = next(pod5_dr.reads(selection=[read_id], preload=['samples']))
        io_read = io.Read.from_pod5_and_alignment(pod5_read, bam_read)


        # Strand filtering
        strand = io_read.ref_reg.strand
        if STRAND_FILTER != "both" and strand != STRAND_FILTER:
            continue

        # --- Q score filtering ---
        q_scores = bam_read.query_qualities
        if np.mean(q_scores) < 22:
            continue

        q_scores_str = clean_array_to_string(q_scores)

        data.append({
            'Read ID': io_read.read_id,
            'Reference Sequence': io_read.ref_seq,
            'Basecalls Length': io_read.seq_len,
            'Reference Mapping Length': io_read.ref_seq_len,
            'Reference Location': str(io_read.ref_reg),
            'Signal': clean_array_to_string(io_read.dacs),
            'Normalized Signal': io_read.norm_signal.tolist(),
            'Ref to Signal': io_read.ref_to_signal.tolist(),
            'Base Sequence': io_read.seq,
            'Sequence to Signal Map': clean_array_to_string(io_read.query_to_signal),
            'Q Scores': q_scores_str,
            'Mapped': io_read.ref_reg is not None

        })

    except Exception as e:
        print(f"⛔ Skipping read {read_id}: {e}")
        continue

df = pd.DataFrame(data)
output_csv = os.path.join(output_dir, 'vis_output.csv')
df.to_csv(output_csv, index=False)
print(f"✅ Saved {len(df)} rows to {output_csv}")

bam_file.close()
pod5_dr.close()


