# xenoquant

xenoquant is a neural network training and inference pipeline for nanopore sequencing of alternative base pairs (XNAs). It builds on Oxford Nanopore Technologies (ONT) workflows and Remora to enable reference-localized detection and quantification of non-standard nucleotides. The toolkit preprocesses FAST5/POD5 data, performs alignment, generates training chunks, and trains recurrent neural network models for XNA-aware basecalling.

xenoquant supports training on defined sequence contexts and performing reference-based XNA identification using trained models. The general workflow consists of two primary steps:

1) Train a model using reference datasets with and without XNA substitutions  
2) Basecall new reads using the trained model  

This software has been tested on ONT R9.4.1 and R10.4.1 flow cells (Flongle and MinION) using Remora 2.x. The project is under active development.

---

## Overview

xenoquant enables:

- Training Remora-compatible RNN models for detecting specific XNA substitutions
- Reference-localized XNA identification
- Export of trained models for use with Dorado, Guppy, or Bonito (via Remora compatibility)
- Handling of non-standard bases using the xFASTA format

Models are trained on specific sequence contexts rather than the full sequence space. For best performance, training should include:

- One dataset containing XNA substitutions
- One dataset containing the corresponding canonical DNA sequence

---

## Command Groups

xenoquant provides two primary commands:

train      Train a xenoquant model on input FAST5 reads using reference sequences containing XNAs.  
basecall   Basecall FAST5 reads around XNA positions using a trained xenoquant model.  

Additional help is available within each command.

---

## Dependencies

xenoquant requires:

- ONT tools (ont-fast5-api, guppy or dorado, Remora 2.x)
- minimap2 (mappy)
- Standard Python scientific packages

A full dependency list is provided in the environment file.

Tested on Ubuntu 18.04 and 20.04.

---

## Training

The `train` command builds a model using two sets of nanopore reads:

1) Reads containing the XNA substitution  
2) Reads containing the canonical DNA comparison  

Each dataset requires a reference FASTA file. Reference FASTA files must contain the XNA base abbreviation in the sequence to define the model focus position.

If both datasets share identical sequence context (with only the XNA position substituted), the same reference file may be reused.

### Training Pipeline

The training workflow follows the Remora 2.x procedure:

1. Convert FAST5 to POD5  
2. Perform initial basecalling  
3. Align reads to reference  
4. Convert FASTA → xFASTA  
5. Generate BED files indicating XNA positions  
6. Generate Remora chunks for XNA and DNA datasets  
7. Merge chunks  
8. Train LTST RNN model (64 base layers)  

Model checkpoints and outputs are written to:

[output_directory]/model

Final model file:

model_best.pt

### Training Command

python xenoquant.py train \
    -w [output_directory] \
    -f [xna_fast5_directory] [dna_fast5_directory] \
    -r [xna_reference.fa] [dna_reference.fa]

---

## Basecalling

The `basecall` command applies a trained model to new reads.

Inputs required:

- Output directory  
- FAST5 directory  
- Reference FASTA (must include XNA positions in sequence)  
- Trained model file (.pt)  

Processing steps:

1. Convert FAST5 → POD5  
2. Perform initial basecalling  
3. Align reads  
4. Generate BED regions from reference  
5. Generate chunks  
6. Apply trained model to XNA positions  

Outputs include:

- Preprocessing files (.pod5, .bam, .bed)  
- Full run summary file  
- Per-read results file (TSV)  

The per-read TSV includes:

- read_id  
- alignment position  
- reference label (1 = XNA, 0 = DNA)  
- predicted class (1 = XNA, 0 = DNA)  
- class probabilities  

### Basecalling Command

python xenoquant.py basecall \
    -w [output_directory] \
    -f [fast5_directory] \
    -r [reference.fa] \
    -m [model_file.pt]

---

## Model Characteristics

xenoquant generates LTST RNN models trained on ±50 signal datapoints around the XNA position defined in the BED file. The corresponding sequence context typically spans ±5–10 bases.

Important notes:

- Models are context-specific and valid only for sequence contexts represented in training data  
- Models do not generalize across flow cell chemistries (e.g., R9.4.1 vs R10.4.1)  
- Validation on new sequence contexts is strongly recommended  
- Models may be converted for use with other ONT basecallers via Remora  

---

## xFASTA Format

Many sequence tools do not support non-ATGC bases. xenoquant uses the xFASTA format to encode XNA positions in the FASTA header.

### Standard FASTA with XNA

>read_header  
ATGGCAACAGGATGABAAGGACGTA

### xFASTA Format

>read_header+X_POS[B:18]  
ATGGCAACAGGATGAGAAGGACGTA

In xFASTA:

- The XNA position is stored in the header  
- The sequence line contains a substituted canonical base  

Default substitutions:

- B → G  
- S → C  
- P → G  
- Z → G  

Substitution rules can be modified in:

lib/xr_params.py

Allowed XNA abbreviations (default):

B, S, P, Z, X, K, J, V

Default pairing rules:

B:S  
P:Z  
X:K  
J:V  

---

## Parameter Customization

Default parameters, paths, and base substitutions can be modified in:

lib/xr_params.py

Users may adjust:

- Confounding base substitutions  
- Model window sizes  
- Allowed XNA bases  
- Workflow defaults  

---

## Example Workflows

### Train a Model

python xenoquant.py train \
    -w path/to/output_directory \
    -f path/to/xna_fast5_directory path/to/dna_fast5_directory \
    -r path/to/xna_reference.fa path/to/dna_reference.fa

Outputs:

- Preprocessed files  
- Model checkpoints  
- model_best.pt in output_directory/model  

Note: Both FASTA references must contain the XNA base in the sequence field.

---

### Basecall Using a Trained Model

python xenoquant.py basecall \
    -w path/to/output_results_directory \
    -f path/to/fast5_directory \
    -r path/to/reference.fa \
    -m path/to/training_directory/model/model_best.pt

Note:

- FASTA reference must contain the XNA base to define evaluation positions  
- Initial basecalling is handled by Guppy or Dorado  

---

## Status

xenoquant is under active development. Use for research and testing purposes.