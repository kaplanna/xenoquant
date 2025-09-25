########################################################################
########################################################################
"""
xr_fasta2x_rc.py 

Title: Unpublished work

Convert FASTA to xFASTA format for XNA base training. 
Handles substitution of modified bases with canonical confounding bases
and generates reverse-complement FASTAs when needed.

By: H. Kawabe, N. Kaplan, J. A. Marchand
Updated: 09/22/25
"""
########################################################################
########################################################################

import sys
import os
from Bio.Seq import Seq
from xr_params import *
from xr_tools import fetch_xna_pos

########################################################################
print("Xemora [Status] - Initializing Xemora...")

# Helpers
########################################################################

def get_confounding_base(base, mod_base, can_base):
    """Return canonical substitute for modified base, or original base if not modified."""
    return can_base if base == mod_base else base

def get_xna_rc(base, xna_pairs):
    """Return reverse complement for an XNA base using the xna_base_pairs list."""
    for pair in xna_pairs:
        if base == pair[0]:
            return pair[1]
        elif base == pair[1]:
            return pair[0]
    # Default DNA complements
    dna_rc = {"A":"T","T":"A","G":"C","C":"G","N":"N"}
    return dna_rc.get(base, base)


########################################################################
# Input / Output
########################################################################

if len(sys.argv) != 3:
    sys.exit("[xFASTA] Usage: python xr_fasta2x_rc.py <input.fasta> <output.fasta>")

input_fasta = sys.argv[1]
output_fasta = sys.argv[2]

detected_xfasta_header = False
detected_xna = False
require_rc_fasta = False

########################################################################
# Conversion
########################################################################

with open(output_fasta, "w") as fout, open(input_fasta, "r") as fin:
    for line in fin:
        if line.startswith(">"):
            header = line.strip()
            if "XPOS[" in header:
                detected_xfasta_header = True
        else:
            seq = line.strip().upper()

            # Detect non-standard bases
            diff = list(set(seq) - set(standard_bases))
            if len(diff) > 0:
                # Handle each XNA type found in sequence
                for xna_base in diff:
                    x_locs = [i for i, c in enumerate(seq) if c == xna_base]

                    # Substitution with canonical base
                    substitution_base = get_confounding_base(xna_base, mod_base, can_base)

                    # Track if this requires RC fasta (segmentation-specific bases)
                    if xna_base in xna_segmentation_model_sets:
                        require_rc_fasta = True

                    # Prepare header
                    header_x = header + "+XPOS[" + "".join(f"{xna_base}:{pos}-" for pos in x_locs)
                    header_x = header_x[:-1] + "]\n"

                    # Prepare substituted sequence
                    clean_seq = seq.replace(xna_base, substitution_base)

                    fout.write(header_x)
                    fout.write(clean_seq + "\n")

                    if write_gaps:
                        gap_header = header + "+-+_GAP[]\n"
                        gap_seq = seq.replace(xna_base, "-").replace("-", "")
                        fout.write(gap_header)
                        fout.write(gap_seq + "\n")

                detected_xna = True

            elif len(diff) == 0 and write_no_xna_seq:
                fout.write(header + "\n")
                fout.write(seq + "\n")

########################################################################
# Reverse complement FASTA (if needed)
########################################################################

if require_rc_fasta:
    if output_fasta.endswith(".fasta"):
        rc_name = output_fasta[:-6] + "_rc.fasta"
    elif output_fasta.endswith(".fa"):
        rc_name = output_fasta[:-3] + "_rc.fa"
    else:
        rc_name = output_fasta + "_rc.fa"

    with open(rc_name, "w") as frc, open(output_fasta, "r") as fo:
        x_pos_to_rc = []
        for line in fo:
            if line.startswith(">") and "GAP" not in line:
                header = line.strip()
                x_pos_base = fetch_xna_pos(header)
                x_pos_to_rc = [[x[0], x[1].replace("]", "")] for x in x_pos_base]

            elif not line.startswith(">") and x_pos_to_rc:
                seq = line.strip()
                seq_rc = str(Seq(seq).reverse_complement())
                seq_len = len(seq_rc)

                header_rc = header.split("+XPOS[")[0] + "+RC+XPOS["

                for base, pos in x_pos_to_rc:
                    rc_base = get_xna_rc(base, xna_base_pairs)
                    rc_sub = get_confounding_base(rc_base, mod_base, can_base)
                    rc_pos = seq_len - int(pos) - 1
                    seq_rc = seq_rc[:rc_pos] + rc_sub + seq_rc[rc_pos+1:]
                    header_rc += f"{rc_base}:{rc_pos-1}-"

                header_rc = header_rc[:-1] + "]\n"
                frc.write(header_rc)
                frc.write(seq_rc + "\n")
else:
    try:
        rmv = output_fasta.replace(".fa", "_rc") + ".fa"
        os.remove(rmv)
    except:
        pass

########################################################################
# Sanity checks
########################################################################

if detected_xfasta_header:
    print("Xenomorph Status - [Error] Fasta input file already in xfasta format")
elif not detected_xna:
    print("Xenomorph Status - [Error] No XNAs detected in fasta input sequence.")
else:
    print("[xFASTA] Conversion complete.")

