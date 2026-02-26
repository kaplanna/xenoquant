# xenoquant

xenoquant is a neural network training and analysis pipeline for nanopore sequencing of alternative base pairs (XNAs). It integrates Oxford Nanopore Technologies (ONT) workflows with Remora-based model training, reference-localized XNA detection, demultiplexing, signal-level visualization, and downstream analysis.

The system consists of two layers:

1) Core engine (`xenoquant.py`) — handles model training and XNA-aware basecalling  
2) Master orchestration script (`xenoquant_pipe.py`) — controls training, basecalling, demultiplexing, analysis, and visualization workflows  

This software has been tested on ONT R9.4.1 and R10.4.1 flow cells (Flongle and MinION) using Remora 2.x. The project is under active development.

---

## Overview

xenoquant supports:

- Training Remora-compatible RNN models for XNA detection
- Reference-localized XNA basecalling
- POD5-based nanopore preprocessing
- Alignment result extraction
- Barcode demultiplexing (cutadapt-based)
- Raw basecall class analysis
- Signal-level visualization (metrics, spaghetti plots, step plots, violin plots)
- Export of trained models for use with Dorado, Guppy, or Bonito (via Remora compatibility)

Models are trained on defined sequence contexts rather than the entire sequence space. For best performance, training should include:

- One dataset containing XNA substitutions
- One dataset containing the corresponding canonical DNA sequence

---

## Repository Structure

Core scripts:

xenoquant.py  
    Command-line interface for model training and basecalling  

xenoquant_pipe.py  
    Master workflow controller with feature switches  

lib/  
    xr_signal_metrics.py  
    xr_signal_plot_v2.py  
    xr_signal_plot_step.py  
    xr_violin.py  
    xr_extract_metrics.py  
    xr_results.py  
    xr_demux.py  
    xr_raw_basecall_analysis.py  
    xr_tools.py  
    xr_params.py  

---

## Dependencies

xenoquant requires:

- ONT tools (guppy or dorado, Remora 2.x)
- minimap2 (mappy)
- cutadapt (for demultiplexing)
- Standard Python scientific libraries

A full dependency list is provided in the environment file.

Tested on Ubuntu 18.04 and 20.04.

---

## Core Commands (xenoquant.py)

### Train

The `train` command builds a model using two nanopore datasets:

1) Reads containing the XNA substitution  
2) Reads containing the canonical DNA comparison  

Each dataset requires a reference FASTA file containing the XNA base abbreviation in the sequence to define the model focus position.

Training workflow:

1. Convert FAST5 (if necessary) to POD5  
2. Perform initial basecalling  
3. Align reads to reference  
4. Convert FASTA → xFASTA  
5. Generate BED files indicating XNA positions  
6. Generate Remora chunks  
7. Merge chunks  
8. Train LTST RNN model (64 base layers)  

Model outputs are written to:

[output_directory]/model/

Final model file:

model_best.pt

Command:

python xenoquant.py train \
    -w [output_directory] \
    -f [xna_pod5_directory] [dna_pod5_directory] \
    -r [xna_reference.fa] [dna_reference.fa]

---

### Basecall

The `basecall` command applies a trained model to new reads.

Inputs:

- Output directory  
- POD5 directory  
- Reference FASTA (must include XNA positions in sequence)  
- Trained model (.pt)  

Processing:

1. Initial basecalling  
2. Alignment  
3. BED generation  
4. Chunk extraction  
5. Model inference  

Outputs:

- Preprocessing files (.pod5, .bam, .bed)  
- Full run summary  
- Per-read TSV results  

Per-read results include:

- read_id  
- alignment position  
- reference label (1 = XNA, 0 = DNA)  
- predicted class (1 = XNA, 0 = DNA)  
- class probabilities  

Command:

python xenoquant.py basecall \
    -w [output_directory] \
    -f [pod5_directory] \
    -r [reference.fa] \
    -m [model_best.pt]

---

## Master Pipeline (xenoquant_pipe.py)

The master pipeline script provides switch-controlled orchestration of all workflow stages.

Feature switches include:

train_model  
basecall_reads  
output_alignment_results  
cutadapt_demux  
raw_basecall_analysis  

Visualization switches:

plot_signal_metrics  
plot_signal_spaghetti  
plot_signal_step  
plot_signal_violin  
extract_metrics  

Each stage is executed via subprocess calls to the appropriate module.

This script allows reproducible execution of full workflows without manually invoking each component.

---

## Demultiplexing

Barcode-based demultiplexing is supported using cutadapt.

Inputs:

- Working directory  
- Barcode pair CSV  

Handled by:

lib/xr_demux.py

Activated via:

cutadapt_demux = True

---

## Raw Basecall Analysis

Performs class-based filtering and summary analysis of basecalling results.

Handled by:

lib/xr_raw_basecall_analysis.py

Activated via:

raw_basecall_analysis = True

---

## Visualization Suite

Signal-level visualizations are supported for trained or basecalled datasets.

Scripts:

xr_signal_metrics.py  
xr_signal_plot_v2.py  
xr_signal_plot_step.py  
xr_violin.py  
xr_extract_metrics.py  

These modules operate on the working directory and require specification of the XNA base (retrieved from xr_params.py).

Visualizations include:

- Signal summary metrics
- Spaghetti plots
- Step-aligned signal plots
- Violin distributions
- Extracted quantitative metrics

---

## Model Characteristics

Models are trained on ±50 signal datapoints surrounding the XNA position defined in the BED file. This corresponds roughly to ±5–10 nucleotide sequence context, depending on pore chemistry.

Important:

- Models are sequence-context specific  
- Models do not generalize across flow cell chemistries  
- Validation on new contexts is recommended  
- Remora compatibility allows conversion to other ONT basecallers  

---

## xFASTA Format

xenoquant uses the xFASTA format to encode XNA positions in FASTA headers.

Standard FASTA with XNA:

>read_header  
ATGGCAACAGGATGABAAGGACGTA  

xFASTA format:

>read_header+X_POS[B:18]  
ATGGCAACAGGATGAGAAGGACGTA  

In xFASTA:

- The header stores XNA position  
- The sequence line contains substituted canonical base  

Default substitutions:

B → A
S → T
P → G  
Z → C

Substitution rules are configurable in:

lib/xr_params.py  

Allowed XNA abbreviations (default):

B, S, P, Z, D, X 

Default pairing rules:

B:S  
P:Z  
Ds:Px

---

## Parameter Customization

Workflow parameters, default paths, base substitutions, and model settings are configurable in:

lib/xr_params.py

---

## Status

xenoquant is under active development and intended for research use.