#!/usr/bin/env bash

# This script prepares paired-end RNA-seq reads for downstream expression analysis.
# The raw FASTQ files contain the sequenced RNA-derived reads together with adapter
# sequences and low-quality bases that can interfere with downstream alignment.
#
# Raw FASTQ files are first collected in a separate input folder to keep the project
# structure organized. The script then checks that every R1 file has a matching R2 file
# before any trimming starts.
#
# fastp is used to trim adapters and perform basic quality filtering while keeping
# the R1 and R2 reads paired. Each sample is processed separately and the cleaned
# FASTQ files are saved for the following RNA-seq alignment and quantification steps.
# fastp also creates HTML and JSON reports that summarize read quality before and
# after trimming. Samples with complete fastp outputs are skipped when the script
# is run again.


# Stops the script if a command fails, an undefined variable is used,
# or a command within a pipeline fails.
set -euo pipefail


# Defines the main project folder and the folders used for raw and trimmed FASTQ files.
ROOT="/home1/RNAseq"

RAW_DIR="${ROOT}/00_raw_fastq"

OUTPUT_DIR="${ROOT}/01_fastp_results"


# Number of processor threads used by fastp to process the reads in parallel.
THREADS=16


# Creates the raw-data and fastp output folders if they do not already exist.
mkdir -p "${RAW_DIR}" "${OUTPUT_DIR}"


# Makes filename patterns with no matches return an empty list instead of
# keeping the unmatched pattern as text.
shopt -s nullglob


# Finds raw R1 and R2 FASTQ files that are still located directly in the project folder.
root_fastq_files=(
    "${ROOT}"/*_R1.fastq.gz
    "${ROOT}"/*_R2.fastq.gz
)


# Moves raw FASTQ files into the dedicated raw-data folder.
# This keeps the original sequencing files separate from analysis outputs.
if (( ${#root_fastq_files[@]} > 0 )); then

    echo "Moving raw FASTQ files to ${RAW_DIR}"


    # Moves each FASTQ file separately and stops if a file with the same
    # name already exists in the destination folder.
    for fastq in "${root_fastq_files[@]}"; do

        destination="${RAW_DIR}/$(basename "${fastq}")"


        # Prevents an existing raw FASTQ file from being overwritten accidentally.
        if [[ -e "${destination}" ]]; then

            echo "ERROR: Destination already exists: ${destination}" >&2
            exit 1
        fi


        mv "${fastq}" "${RAW_DIR}/"

    done

else

    # Uses the existing raw-data folder if the FASTQ files were already moved earlier.
    echo "No raw FASTQ files found directly in ${ROOT}; using ${RAW_DIR}."

fi


# Finds all R1 FASTQ files, which are used to identify the paired-end samples.
r1_files=("${RAW_DIR}"/*_R1.fastq.gz)


# Stops the script if no R1 FASTQ files matching the expected naming pattern were found.
if (( ${#r1_files[@]} == 0 )); then

    echo "ERROR: No *_R1.fastq.gz files found in ${RAW_DIR}" >&2
    exit 1

fi


# Checks all samples before trimming starts to make sure every R1 file
# has a corresponding R2 file.
for R1 in "${r1_files[@]}"; do

    # Extracts the sample name from the R1 filename.
    SAMPLE=$(basename "${R1}" _R1.fastq.gz)


    # Uses the sample name from R1 to define the corresponding R2 file.
    R2="${RAW_DIR}/${SAMPLE}_R2.fastq.gz"


    # Stops the script if the matching R2 file is missing or empty because
    # paired-end reads should be processed together.
    [[ -s "${R2}" ]] || {
        echo "ERROR: Missing R2 for ${SAMPLE}: ${R2}" >&2
        exit 1
    }

done


# Processes each paired-end RNA-seq sample separately.
for R1 in "${r1_files[@]}"; do

    # Extracts the sample name from the R1 filename.
    SAMPLE=$(basename "${R1}" _R1.fastq.gz)


    # Uses the sample name from R1 to define the corresponding R2 file.
    R2="${RAW_DIR}/${SAMPLE}_R2.fastq.gz"


    # Defines the trimmed FASTQ files and fastp quality-control reports
    # that are created for the current sample.
    OUT_R1="${OUTPUT_DIR}/${SAMPLE}_trimmed_R1.fastq.gz"
    OUT_R2="${OUTPUT_DIR}/${SAMPLE}_trimmed_R2.fastq.gz"
    HTML="${OUTPUT_DIR}/${SAMPLE}_fastp.html"
    JSON="${OUTPUT_DIR}/${SAMPLE}_fastp.json"


    # Skips the sample if both trimmed FASTQ files and both quality-control
    # reports already exist and are not empty.
    if [[ -s "${OUT_R1}" && -s "${OUT_R2}" && -s "${HTML}" && -s "${JSON}" ]]; then

        echo "Skipping ${SAMPLE}: complete fastp outputs already exist."
        continue

    fi


    echo "Processing sample: ${SAMPLE}"


    # Runs fastp on the paired-end reads. Adapter sequences are detected automatically,
    # and the cleaned R1 and R2 reads are written to separate trimmed FASTQ files.
    # HTML and JSON reports store quality-control information for each sample.
    fastp \
        --in1 "${R1}" \
        --in2 "${R2}" \
        --out1 "${OUT_R1}" \
        --out2 "${OUT_R2}" \
        --thread "${THREADS}" \
        --detect_adapter_for_pe \
        --html "${HTML}" \
        --json "${JSON}"


    echo "Finished: ${SAMPLE}"

done


# Confirms that all required fastp processing is finished.
echo "All fastp processing finished."