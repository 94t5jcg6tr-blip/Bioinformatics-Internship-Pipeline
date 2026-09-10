#!/bin/bash

# This script aligns trimmed paired-end reads to a combined mm10/dm6 reference genome.
# The Drosophila reads come from the spike-in control and are kept separately later
# so they can be used for normalization between ChIP-seq samples.
# Bowtie2 is used for the alignment, followed by Samtools to prepare the BAM files
# for duplicate removal.
#
# The reads are first sorted by name and processed with fixmate to add information
# about the paired read. They are then sorted by genomic position and split into
# mouse (mm10) and Drosophila (dm6) reads. Duplicates are removed separately from
# both BAM files and the final files are indexed.


# Stop the script on errors and unset variables.
set -Eeuo pipefail


# Make sure the correct Conda environment is active.
if [[ "${CONDA_DEFAULT_ENV:-}" != "NGS-py37" ]]; then
    echo "ERROR: Please activate NGS-py37 before running this script."
    exit 1
fi


# Bowtie2 uses more threads for the alignment step.
THREADS=60

# Samtools uses fewer threads for sorting and BAM processing.
SORT_THREADS=8


# Combined mouse and Drosophila Bowtie2 reference.
reference_file="/home1/Bowtie2_Index_mm10withdm6/mm10withdm6"
reference_name=$(basename "$reference_file")


# Check that Bowtie2 is available.
if ! command -v bowtie2 >/dev/null 2>&1; then
    echo "Error: bowtie2 not found in NGS-py37."
    exit 1
fi


# Check that Samtools is available.
if ! command -v samtools >/dev/null 2>&1; then
    echo "Error: samtools not found in NGS-py37."
    exit 1
fi


# Check that the Bowtie2 reference index exists.
if [[ ! -f "${reference_file}.1.bt2" && ! -f "${reference_file}.1.bt2l" ]]; then
    echo "Error: Bowtie2 index not found: ${reference_file}"
    exit 1
fi


echo "Running Bowtie2/Samtools in environment: $CONDA_DEFAULT_ENV"
echo "Reference genome: $reference_file"
echo "Bowtie2 threads: $THREADS"
echo "Samtools sort threads: $SORT_THREADS"


# Create folders for BAM files, logs, QC and temporary files.
mkdir -p 02_alignment/bam_combined
mkdir -p 02_alignment/bam_mm10
mkdir -p 02_alignment/bam_dm6
mkdir -p 02_alignment/logs
mkdir -p 02_alignment/qc
mkdir -p 02_alignment/tmp


# Get sample names from the trimmed R1 FASTQ files.
sample_names=$(ls 01_fastp/trimmed_fastq/*_trimmed_R1.fastq.gz 2>/dev/null | \
    sed -e 's|01_fastp/trimmed_fastq/||g' -e 's/_trimmed_R1.fastq.gz//g' || true)


# Stop if no trimmed samples were found.
if [ -z "$sample_names" ]; then
    echo "Error: No trimmed R1 files found in 01_fastp/trimmed_fastq/"
    exit 1
fi


echo "Samples found:"
echo "$sample_names"


# Process each sample separately.
for input_file in $sample_names
do

    R1="01_fastp/trimmed_fastq/${input_file}_trimmed_R1.fastq.gz"
    R2="01_fastp/trimmed_fastq/${input_file}_trimmed_R2.fastq.gz"

    # Set output and log files for this sample.
    COMBINED_BAM="02_alignment/bam_combined/${input_file}_${reference_name}_sorted.bam"
    MM10_BAM="02_alignment/bam_mm10/${input_file}_mm10_sorted_markdup.bam"
    DM6_BAM="02_alignment/bam_dm6/${input_file}_dm6_sorted_markdup.bam"
    LOG="02_alignment/logs/Bowtie2_log_${input_file}.txt"


    # Skip samples if both final BAM files already exist.
    if [[ -s "$MM10_BAM" && -s "$DM6_BAM" ]]; then
        echo "Skipping ${input_file}: final BAM files already exist."
        continue
    fi


    # Check that both paired-end FASTQ files are available.
    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
        echo "Error: Missing trimmed R1 or R2 for ${input_file}"
        exit 1
    fi


    echo "Running Bowtie2/Samtools for ${input_file}..."


    # Align paired-end reads to the combined mm10/dm6 reference.
    bowtie2 -p "$THREADS" -x "$reference_file" \
        --no-mixed --no-discordant --no-unal \
        -1 "$R1" \
        -2 "$R2" \
        2> >(tee "$LOG" >&2) | \

    # Sort reads by name so paired reads are next to each other for fixmate.
    samtools sort -@ "$SORT_THREADS" -n \
        -T "02_alignment/tmp/${input_file}_namesort" - | \

    # Add mate information, such as the position of the paired read,
    # which is needed for duplicate detection.
    samtools fixmate -m - - | \

    # Sort reads by genomic position before saving the combined BAM file.
    samtools sort -@ "$SORT_THREADS" \
        -T "02_alignment/tmp/${input_file}_positionsort" \
        -o "$COMBINED_BAM"


    # Index the combined coordinate-sorted BAM file.
    samtools index "$COMBINED_BAM"


    echo "Splitting ${input_file} into mm10/dm6 and removing duplicates..."


    # Extract mouse reads for the main ChIP-seq analysis.
    samtools view -@ "$SORT_THREADS" -b "$COMBINED_BAM" \
        $(printf 'chr%s ' {1..19}) chrX chrY | \

    # Sort mouse reads again before duplicate removal.
    samtools sort -@ "$SORT_THREADS" \
        -T "02_alignment/tmp/${input_file}_mm10sort" - | \

    # Remove duplicate mouse reads.
    samtools markdup -@ "$SORT_THREADS" -r - "$MM10_BAM"


    # Extract Drosophila spike-in reads for later normalization.
    samtools view -@ "$SORT_THREADS" -b "$COMBINED_BAM" \
        dm6_chr2L dm6_chr2R dm6_chr3L dm6_chr3R dm6_chr4 \
        dm6_chrX dm6_chrY dm6_chrM | \

    # Sort Drosophila reads before duplicate removal.
    samtools sort -@ "$SORT_THREADS" \
        -T "02_alignment/tmp/${input_file}_dm6sort" - | \

    # Remove duplicate Drosophila reads.
    samtools markdup -@ "$SORT_THREADS" -r - "$DM6_BAM"


    # Index the final mm10 and dm6 BAM files.
    samtools index "$MM10_BAM"
    samtools index "$DM6_BAM"


    # Save read and mapping counts for the mm10 and dm6 BAM files.
    samtools flagstat "$MM10_BAM" \
        > "02_alignment/qc/${input_file}_mm10_flagstat.txt"

    samtools flagstat "$DM6_BAM" \
        > "02_alignment/qc/${input_file}_dm6_flagstat.txt"


    echo "Bowtie2 and Samtools finished for ${input_file}."

done


echo "All Bowtie2/Samtools jobs finished."