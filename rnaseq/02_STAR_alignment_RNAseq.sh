#!/usr/bin/env bash

# This script aligns the trimmed paired-end RNA-seq reads to the mouse reference genome.
# After adapter trimming, the reads still need to be assigned to their genomic positions
# before gene expression can be quantified and compared between the different genotypes.
#
# STAR is used as the RNA-seq aligner because it can map reads across exon-exon junctions,
# which is important for transcript-derived sequencing reads. Each sample is processed
# separately using the same pre-built STAR genome index.
#
# The aligned reads are saved as coordinate-sorted BAM files, and STAR also generates
# gene-level read counts that can be used for downstream expression analysis.
# Samples with complete STAR outputs are skipped when the script is run again, and
# the final BAM files are indexed with samtools for faster downstream access.


# Stops the script if a command fails, an undefined variable is used,
# or a command within a pipeline fails.
set -euo pipefail


# Defines the main project folder and the input, output and reference locations.
ROOT="/home1/RNAseq"

TRIMMED_DIR="${ROOT}/01_fastp_results"

OUT_DIR="${ROOT}/02_alignment"


# Folder containing the pre-built STAR index of the reference genome.
# The index allows STAR to map reads to the genome efficiently.
STAR_INDEX="/home1/RSEM-STAR_Index_v2-7-10b"


# Number of processor threads used by STAR and samtools.
THREADS=60


# Creates the alignment output folder if it does not already exist.
mkdir -p "${OUT_DIR}"


# Stops the script if the trimmed-read folder is missing.
[[ -d "${TRIMMED_DIR}" ]] || {
    echo "ERROR: Missing ${TRIMMED_DIR}" >&2
    exit 1
}


# Stops the script if the STAR reference index folder is missing.
[[ -d "${STAR_INDEX}" ]] || {
    echo "ERROR: Missing STAR index ${STAR_INDEX}" >&2
    exit 1
}


# Makes filename patterns with no matches return an empty list instead of
# keeping the unmatched pattern as text.
shopt -s nullglob


# Finds all trimmed R1 FASTQ files, which are used to identify the RNA-seq samples.
r1_files=("${TRIMMED_DIR}"/*_trimmed_R1.fastq.gz)


# Stops the script if no trimmed R1 FASTQ files were found.
if (( ${#r1_files[@]} == 0 )); then

    echo "ERROR: No trimmed R1 FASTQ files found in ${TRIMMED_DIR}" >&2
    exit 1

fi


# Checks all samples before alignment starts to make sure every R1 file
# has a corresponding non-empty R2 file.
for R1 in "${r1_files[@]}"; do

    # Extracts the sample name from the trimmed R1 filename.
    SAMPLE=$(basename "${R1}" _trimmed_R1.fastq.gz)


    # Uses the sample name from R1 to define the corresponding trimmed R2 file.
    R2="${TRIMMED_DIR}/${SAMPLE}_trimmed_R2.fastq.gz"


    # Stops the script if the matching R2 file is missing or empty because
    # paired-end reads should be aligned together.
    [[ -s "${R2}" ]] || {
        echo "ERROR: Missing trimmed R2 for ${SAMPLE}: ${R2}" >&2
        exit 1
    }

done


# Processes each trimmed paired-end RNA-seq sample separately.
for R1 in "${r1_files[@]}"; do

    # Extracts the sample name from the trimmed R1 filename.
    SAMPLE=$(basename "${R1}" _trimmed_R1.fastq.gz)


    # Uses the sample name from R1 to define the corresponding trimmed R2 file.
    R2="${TRIMMED_DIR}/${SAMPLE}_trimmed_R2.fastq.gz"


    # Defines the common STAR output prefix for the current sample.
    PREFIX="${OUT_DIR}/${SAMPLE}_"


    # Defines the main STAR output files used to check whether the sample
    # has already been processed completely.
    BAM="${PREFIX}Aligned.sortedByCoord.out.bam"

    FINAL_LOG="${PREFIX}Log.final.out"


    # Skips the sample if the coordinate-sorted BAM file and final STAR log
    # already exist and are not empty.
    if [[ -s "${BAM}" && -s "${FINAL_LOG}" ]]; then

        echo "Skipping ${SAMPLE}: complete STAR outputs already exist."
        continue

    fi


    echo "Processing sample: ${SAMPLE}"


    # Runs STAR to align the paired-end reads to the reference genome.
    # Compressed FASTQ files are decompressed during reading, the aligned reads
    # are written directly as a coordinate-sorted BAM file, and STAR also
    # generates gene-level read counts for downstream expression analysis.
    STAR \
        --runThreadN "${THREADS}" \
        --genomeDir "${STAR_INDEX}" \
        --readFilesIn "${R1}" "${R2}" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${PREFIX}" \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts


    # Creates an index for the sorted BAM file so genomic regions can be
    # accessed directly without reading through the complete BAM file.
    samtools index \
        -@ "${THREADS}" \
        "${BAM}"


    echo "Finished: ${SAMPLE}"

done


# Confirms that all required STAR alignments are finished.
echo "All STAR alignments finished."