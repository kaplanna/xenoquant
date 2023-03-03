```                                                                                  

██╗  ██╗███████╗███╗   ███╗ ██████╗ ██████╗  █████╗ 
╚██╗██╔╝██╔════╝████╗ ████║██╔═══██╗██╔══██╗██╔══██╗
 ╚███╔╝ █████╗  ██╔████╔██║██║   ██║██████╔╝███████║
 ██╔██╗ ██╔══╝  ██║╚██╔╝██║██║   ██║██╔══██╗██╔══██║
██╔╝ ██╗███████╗██║ ╚═╝ ██║╚██████╔╝██║  ██║██║  ██║
╚═╝  ╚═╝╚══════╝╚═╝     ╚═╝ ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝
                                                    
```

# Xemora : An XNA sequencing neural network trainer 

## About 
Xemora is a tool pipeline used for nanopore sequencing of alternative basepairs (XNAs) that latches onto Remora 2.0 (ONT). This toolkit incorporates ONT-workflows to preprocess fast5, pod5, fastq, bam, and bed files for training a remora model. XNA sequence handling is done using the xFASTA format. Therefore, FASTA reference inputs should contain XNA bases in sequence lines. Models for basecalling can be trained on specific set of sequences (rather than entire sequence space). For optimal implementation, xemora should be trained on at least two datasets, with and without the XNA substitutions. Xemora models can be exported and used as remora models for guppy, dorado, or bonito basecalling. Alternatively, xemora can handle reference-based XNA identification directly. The general pipeline consists of two steps: 1) Training xemora on a reference set of reads with and without XNAs 2) basecalling using the trained model. 

Xemora command groups (additional help available within each command group):
	train		[Train] a xemora model on a set of input fast5 reads, localized to reference fasta containing XNAs. 
	basecall	[Basecall] a fast5 reads around XNA using a previously trained xemora model. 


This version of Xenomorph was built and tested on Oxford Nanopore Technologies r9.4.1 and r10.4.1 flow cells (Flongle/MinION) using Remora 2.0. xemora is currently under development and for testing purposes only.

This repository is maintained by the XenoBiology Research Group at the University of Washington, Department of Chemical Engineering. 

## Sample nanopore sequences with XNA basepairs 
This toolkit was created to work with a series of non-standard nucleotides that can form the basis of an expanded genetic alphabet (up to 12 letters). In addition to the standard base pairs (A:T, G:C), the XNA models described in this work allow for single xenonucleotide detection of specific forms of B:S, P:Z, X:K, and J:V. Sample nanopore XNA sequences with each of these base pairs encoded can be downloaded from the Sequence Reads Archive (SRA Bioproject: PRJNA932328). More information can be found with the associate publication.

## Dependencies
Xemora requires ONT tools (ont-fast5-api, tombo, guppy, remora 2.0), minimap2 (mappy), and various python packages. A full list of dependencies can be found in the xemora-env.yml document. To use conda for installing dependencies, simply load xemora-env.yml into a new environment. Xemora was built and tested on Ubuntu 18.04 and 20.04. 

        conda env create -f xemora-env.yml

To enter xemora conda environment, then use: 

        conda activate xemora-env


## Xemora command groups 
	train		[Train] a xemora model on a set of input fast5 reads, localized to reference fasta containing XNAs. 
	basecall	[Basecall] a fast5 reads around XNA using a previously trained xemora model. 


### Train
Xemora train takes raw nanopore multi-FAST5 and reference sequences in FASTA format as inputs. To train a xemora model, two sets of fast5 reads/fasta references are required. 1) a set of fast5 reads that contain an XNA with a reference fasta file contains the location of the XNA (indicated by XNA base abbreviation BSPZXKJV); 2) a set of fast5 reads that contain a standard nucleic acid with a reference fasta file that contains the location of comparison (indicated by XNA base abbreviation BSPZXKJV). Though two sets of fasta files are required, model training can be performed from data in a mixed run if reads are properly barcoded. In this event, the same fast5 file can be used as the input. If the two set of data contain the same sequences (but with the XNA position substituted), the same reference file can be used as input. Note: Reference fasta for standard DNA sequences requires XNA sequences in the location of training for specifying model focus region. 

The full training pipeline follows Remora 2.0 procedure. Fast5 files are first converted to Pod5 format. Pod5 files are basecalled using guppy and reference aliged. Bed files are generated from the reference files following fasta to xfasta conversion to indicate positions for model training. Remora chunks are then generated for the XNA sequences and the standard DNA sequences. Chunks are merged in preparation for training. Training proceeds using LTST RNN model with 64 base layers. Model output with checkpoints can be found in the output directory /model folder. 


        python xemora.py train -w [output directory] -f [xna_fast5_diretory] [dna_fast5_directory] -r [xna reference_fasta] [dna_reference_fasta] 

### Basecall
Xemora basecall uses the remora validate.py to basecall on a set of chunks data. As inputs, users need to specify an output directory, a raw fast5 directory, a reference fasta file (with XNAs in positions to basecall), and a model file (model_file.pt) generated by xemora train. An initial set of basecalls is handled by ONT guppy, which is then used for alignment. The output directory is used to store preprocessing files (.pod5, .bam, .bed) and generate two summary results file: a full run summary file and a per-read results file. The per-read results are stored as a tab seperated file and indicate read id, position in alignment where XNA was found, reference base (1 = XNA, 0 = DNA), basecall (1 = XNA, 0 = DNA), and matching probability for DNA,XNA base at that position.

        python xemora.py basecall -w [output directory] -f [xna_fast5_diretory] -r [xna reference_fasta] -m [model_file.pt]



### X(R)emora Models
Xemora generates LTST RNN models that can either be used by various basecallers. Since xemora is built around reference-based basecalling of XNAs, models generated are not intended to be universallly applied to any sequence context. Chunks used for training are set to +/-50 datapoints before and after the specified XNA position in the reference bed file. The sequence context this corresponds to is variable, but can range from +/- 5-10 bases. Though xemora is trained on a specific sequence context, it does not use standard 'motif driven' model workflow which specifies region of testing. Instead, xemora modifications are treated as a global modification on a specified standard base (e.g., any G could be a P), with the location of testing specified by regions denoted as modified in the reference .bed file. This means xemora models can be trained on any arbitrary set of sequence contexts, modifications, and alternative substitutions but will only be valid for the subset of sequence contexts it was trained on. It is recommended that model basecalling be verified on new sequences where the modification being tested is in a slightly different global context (i.e., changing bases distal to XNA positions). Note that model validity does not extend across flow cell types (r9.4.1 vs r10.4.1), so training data should be ran on the flow cells that will ultimately be used for base calling. Though not recommended since xemora trains on pseudo-motifs, models can be converted to bonito or dorado models for use with other basecallers using remora. 

### About xFASTA format 
Many tools used to read and manipulate nucleic acid sequences are not built to handle non ATGC bases. In an effort to streamline how we handle XNAs, xemora was built to handle non-standard bases in FASTA sequences (e.g. BSPZJVKX). The xFASTA file format (.fa) stores XNA positional information in the header of each sequence. xFASTA files can be generated from standard Fasta files that contain non-ATGC bases (e.g. BSPZJVKX) in the sequence. xFASTA files are automatically generated in the standard xemora preprocessing workflow. The fasta2x command is provided for utility, but generally not required. Note that XNA bases in the input sequence will be replaced for a standard ATGC. Signal matching to sequence is highly sensitive what base is chosen as the replacement base. As default, XNA bases are replaced as followed: B>G, S>C, P>G, Z>G. Base substitution settings can be modified in lib/xr_params.py by changing the paired base in the confounding_base variable. 


        Standard FASTA format with XNA in sequence
                >header_for_a_read
                ATGGCAACAGGATGABAAGGACGTA

        xFASTA format with XNA information stored in header. (B replaced with G in sequence)
                >header_for_a_read+X_POS[B:18]
                ATGGCAACAGGATGAGAAGGACGTA

### Modifying parameters and default paths
Parameters, defaults, and file paths can be tuned or modified for desired application. Allowed XNA base abbreviations are by default B,S,P,Z,X,K,J,V with the following base pairing rules (B:S, P:Z, X:K, J:V). Though default parameters and bases are recommended for most applications, parameters can be modified by editing lib/xr_params.py. 



## Examples

#### Train a xemora model


        1) Preprocess two sets of raw nanopore runs and train model (generates: pod5,fastq, bam, bed, chunks for each then trains model)
            python xemora.py train -w path/to/output_directory -f path/to/xna_fast5_directory path/to/dna_fast5_directory  -r path/to/xna_fasta_reference.fa path/to/dna_fasta_reference.fa 
        *) Xemora trained model (from remora output) will be output to path/to/output_directory/model. 
        **) Model for xemora basecalling will be generated as path/to/output_directory/model/best_model.pt. 
        ***) Both fasta file inputs need to contain XNA of interest in sequence field for specifying location of model training 


#### Basecall using xemora model


        1) Basecall set of reads using specified model (generates: pod5,fastq, bam, bed, chunks then applies trained model around specified locations in bed file)
            python xemora.py basecall -w path/to/output_results_directory -f path/to/fast5_directory -r path/to/fasta_reference -m path/to/training_directory/model/best_model.pt. 
        *) Fasta file input needs to contain XNA of interest in sequence field for specifying location of basecalling. 
        **) Initial basecalling handled by guppy.




## Cite us or read more about this work 
    Title: Unpublished work

    By: H. Kawabe, N. Kaplan, J. A. Marchand

