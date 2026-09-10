#!/bin/bash

# This script quantifies gene and transcript expression from the trimmed paired-end
# RNA-seq reads for all samples, including both 3S1G-OHT replicates.
#
# RSEM estimates how many reads originate from each gene and transcript, which provides
# the expression values used later for differential expression analysis.
# For each sample, the script finds the matching R1 and R2 FASTQ files and uses the
# pre-built RSEM reference for quantification.
#
# STAR is used by RSEM for the alignment step, and the resulting gene-level and
# transcript-level expression estimates are saved separately for every sample.
# Samples with complete RSEM outputs are skipped when the script is run again.


# Stops the script immediately if a command fails.
set -e

# Stops the script if an undefined variable is used.
set -u


# ============================================================
# RNA-seq RSEM quantification for all samples
# Includes 3S1G-OHT replicate 1 and replicate 2
# ============================================================


# -----------------------
# Settings
# -----------------------

# Number of processor threads used by RSEM and STAR during quantification.
THREADS=60


# Defines the main RNA-seq project folder.
PROJECT_DIR="/home1/RNAseq"


# Folder containing the trimmed paired-end FASTQ files.
TRIMMED_DIR="$PROJECT_DIR/01_fastp_results"

# Folder containing the pre-built RSEM reference files.
RSEM_REF_DIR="/home1/RSEM-STAR_Index_v2-7-10b"

# Folder where the RSEM quantification results are saved.
OUT_DIR="$PROJECT_DIR/03_RSEM_results"


# Creates the output folder if it does not already exist.
mkdir -p "$OUT_DIR"


# -----------------------
# Find RSEM reference prefix
# -----------------------

# Makes filename patterns with no matches return an empty list instead of
# keeping the unmatched pattern as text.
shopt -s nullglob


# Finds the .grp file belonging to the RSEM reference.
# This file is part of the reference generated for RSEM and is used here
# to identify the common reference name required by rsem-calculate-expression.
RSEM_GRP_FILES=("$RSEM_REF_DIR"/*.grp)


# Stops the script if no RSEM reference file was found.
if [ ${#RSEM_GRP_FILES[@]} -eq 0 ]; then

    echo "Error: no RSEM .grp reference file found."
    echo "Checked folder: $RSEM_REF_DIR"

    exit 1

fi


# Uses the first .grp file to determine the RSEM reference prefix.
# RSEM expects the shared filename prefix of the reference files rather than
# the .grp file itself, so the .grp extension is removed.
RSEM_REF="${RSEM_GRP_FILES[0]}"

RSEM_REF=${RSEM_REF%.grp}


# Prints the RSEM reference that will be used for quantification.
echo "Using RSEM reference:"
echo "$RSEM_REF"
echo


# -----------------------
# Find trimmed FASTQ files
# -----------------------

# Finds all trimmed R1 FASTQ files, which are used to identify the RNA-seq samples.
R1_FILES=("$TRIMMED_DIR"/*_trimmed_R1.fastq.gz)


# Stops the script if no trimmed R1 files matching the expected pattern were found.
# The folder content is printed to make path or filename problems easier to identify.
if [ ${#R1_FILES[@]} -eq 0 ]; then

    echo "Error: no trimmed R1 FASTQ files found."
    echo "Checked folder: $TRIMMED_DIR"
    echo "Expected pattern: *_trimmed_R1.fastq.gz"
    echo
    echo "Files in this folder are:"

    ls -lh "$TRIMMED_DIR" | head -50

    exit 1

fi


# -----------------------
# Run RSEM
# -----------------------

# Processes each trimmed paired-end RNA-seq sample separately.
for R1 in "${R1_FILES[@]}"

do

    # Extracts the sample name from the trimmed R1 filename.
    SAMPLE=$(basename "$R1" _trimmed_R1.fastq.gz)


    # Uses the sample name from R1 to define the corresponding trimmed R2 file.
    R2="$TRIMMED_DIR/${SAMPLE}_trimmed_R2.fastq.gz"


    # Stops the analysis if the matching R2 file is missing or empty because
    # paired-end reads should be quantified together.
    if [ ! -s "$R2" ]; then

        echo "Error: missing or empty trimmed R2 file for sample $SAMPLE"
        echo "Expected: $R2"

        exit 1

    fi


    # Defines the main gene-level and transcript-level RSEM output files
    # used to check whether this sample has already been quantified completely.
    GENES_RESULT="$OUT_DIR/${SAMPLE}.genes.results"

    ISOFORMS_RESULT="$OUT_DIR/${SAMPLE}.isoforms.results"


    # Skips the sample if both main RSEM result files already exist and are not empty.
    if [ -s "$GENES_RESULT" ] && [ -s "$ISOFORMS_RESULT" ]; then

        echo "Skipping $SAMPLE: complete RSEM outputs already exist."
        continue

    fi


    echo "Processing sample: $SAMPLE"


    # Runs RSEM to estimate gene and transcript expression from the paired-end reads.
    # STAR is used internally for alignment, compressed FASTQ files are read directly,
    # and the resulting expression estimates are saved with the sample-specific prefix.
    rsem-calculate-expression \
        --paired-end \
        --star \
        --star-gzipped-read-file \
        -p "$THREADS" \
        "$R1" "$R2" \
        "$RSEM_REF" \
        "$OUT_DIR/$SAMPLE"


    echo "Finished: $SAMPLE"
    echo

done


# Confirms that all required RSEM quantifications are finished.
echo "All RSEM quantifications finished."