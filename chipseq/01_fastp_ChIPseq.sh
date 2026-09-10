#!/bin/bash

# This script trims paired-end FASTQ files before downstream analysis.
# fastp is used as the main trimming tool. If a sample cannot be processed
# successfully, it is moved to a retry folder and processed again with
# fixed Illumina adapter sequences. Samples that still fail are finally
# processed with cutadapt as a fallback.
#
# Time limits are used so that single problematic FASTQ files do not block
# the complete run. Successful samples are moved to the raw FASTQ archive,
# while failed samples are passed to the next retry step. Logs, reports and
# status files are saved separately for each run.


# Stop the script on errors and unset variables.
set -Eeuo pipefail

# Ignore file patterns if no matching files are found.
shopt -s nullglob


# Number of threads used for fastp and cutadapt.
THREADS=60

# Time limits prevent single problematic files from blocking the whole run.
TIMEOUT_FASTP_NORMAL="45m"
TIMEOUT_FASTP_EXPLICIT="150m"
TIMEOUT_CUTADAPT="300m"


# Illumina TruSeq adapter sequences used for the explicit fallback.
ADAPTER_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"


# Add a timestamp to keep log and status files from different runs separate.
RUN_ID=$(date +%Y%m%d_%H%M%S)


# Check that the required programs are available.
command -v fastp >/dev/null || {
    echo "Error: fastp not found in NGS-py37"
    exit 1
}

command -v timeout >/dev/null || {
    echo "Error: timeout not found"
    exit 1
}


# Install cutadapt if it is needed for the fallback and not available yet.
if ! command -v cutadapt >/dev/null 2>&1; then
    echo "cutadapt not found. Installing..."
    conda install -c bioconda -c conda-forge cutadapt -y
fi


# Create folders for outputs, logs and failed samples.
mkdir -p 01_fastp/trimmed_fastq
mkdir -p 01_fastp/logs
mkdir -p 01_fastp/reports
mkdir -p 01_fastp/status
mkdir -p 00_raw_fastq
mkdir -p 00_retry_fastp_explicit
mkdir -p 00_retry_cutadapt
mkdir -p 00_final_failed


# Set run-specific log and status files.
MAIN_LOG="01_fastp/logs/fastp_cutadapt_run_${RUN_ID}.log"
DONE="01_fastp/status/completed_samples_${RUN_ID}.tsv"
FAILED="01_fastp/status/final_failed_samples_${RUN_ID}.tsv"


# Add headers to the status files.
echo -e "sample\tmethod\tR1_out\tR2_out" > "$DONE"
echo -e "sample\treason" > "$FAILED"


# Write terminal output to the main log file as well.
exec > >(tee -a "$MAIN_LOG") 2>&1


echo "======================================"
echo "Run started: $(date)"
echo "Environment: $CONDA_DEFAULT_ENV"
echo "Threads: $THREADS"
echo "FastP normal timeout: $TIMEOUT_FASTP_NORMAL"
echo "FastP explicit-adapter timeout: $TIMEOUT_FASTP_EXPLICIT"
echo "Cutadapt timeout: $TIMEOUT_CUTADAPT"
echo "======================================"

fastp --version
cutadapt --version



# First trimming attempt with fastp and automatic PE adapter detection.
run_fastp_normal () {

    sample="$1"
    input_dir="$2"
    fail_dir="$3"

    # Set input, output and report files for this sample.
    R1="${input_dir}/${sample}_R1.fastq.gz"
    R2="${input_dir}/${sample}_R2.fastq.gz"

    OUT1="01_fastp/trimmed_fastq/${sample}_trimmed_R1.fastq.gz"
    OUT2="01_fastp/trimmed_fastq/${sample}_trimmed_R2.fastq.gz"

    LOG="01_fastp/logs/FastpLog_normal_${sample}.txt"
    HTML="01_fastp/reports/FastpLog_normal_${sample}.html"
    JSON="01_fastp/reports/FastpLog_normal_${sample}.json"

    echo "--------------------------------------"
    echo "FastP normal: $sample"
    echo "Input dir: $input_dir"
    echo "Started: $(date)"


    # Skip samples that were already processed.
    if [[ -s "$OUT1" && -s "$OUT2" ]]; then
        echo "Skipping $sample: trimmed files already exist."
        return 0
    fi


    # Stop this attempt if one of the paired FASTQ files is missing or empty.
    if [[ ! -s "$R1" || ! -s "$R2" ]]; then
        echo "Missing/empty raw FASTQ for $sample"
        return 1
    fi


    # Remove incomplete output files before rerunning the sample.
    rm -f "$OUT1" "$OUT2" "$HTML" "$JSON"


    # Run fastp with a time limit and save its exit status.
    set +e

    timeout "$TIMEOUT_FASTP_NORMAL" fastp \
        -w "$THREADS" \
        -i "$R1" \
        -I "$R2" \
        -o "$OUT1" \
        -O "$OUT2" \
        --detect_adapter_for_pe \
        -h "$HTML" \
        -j "$JSON" \
        > "$LOG" 2>&1

    status=$?

    set -e


    # Move problematic samples to the next retry step.
    if [[ "$status" -ne 0 || ! -s "$OUT1" || ! -s "$OUT2" ]]; then

        echo "FastP normal failed/timed out for $sample"

        rm -f "$OUT1" "$OUT2"

        mkdir -p "$fail_dir"
        mv "$R1" "$R2" "$fail_dir"/

        return 1
    fi


    # Archive the raw FASTQ files and record the successful run.
    mv "$R1" "$R2" 00_raw_fastq/

    echo -e "${sample}\tfastp_normal\t${OUT1}\t${OUT2}" >> "$DONE"

    echo "Finished FastP normal: $sample"

    return 0
}



# Retry problematic samples using fixed Illumina adapter sequences.
run_fastp_explicit () {

    sample="$1"
    input_dir="$2"
    fail_dir="$3"

    # Set input, output and report files for this sample.
    R1="${input_dir}/${sample}_R1.fastq.gz"
    R2="${input_dir}/${sample}_R2.fastq.gz"

    OUT1="01_fastp/trimmed_fastq/${sample}_trimmed_R1.fastq.gz"
    OUT2="01_fastp/trimmed_fastq/${sample}_trimmed_R2.fastq.gz"

    LOG="01_fastp/logs/FastpLog_explicit_${sample}.txt"
    HTML="01_fastp/reports/FastpLog_explicit_${sample}.html"
    JSON="01_fastp/reports/FastpLog_explicit_${sample}.json"

    echo "--------------------------------------"
    echo "FastP explicit adapters: $sample"
    echo "Input dir: $input_dir"
    echo "Started: $(date)"


    # Skip samples that were already processed.
    if [[ -s "$OUT1" && -s "$OUT2" ]]; then
        echo "Skipping $sample: trimmed files already exist."
        return 0
    fi


    # Stop this attempt if one of the paired FASTQ files is missing or empty.
    if [[ ! -s "$R1" || ! -s "$R2" ]]; then
        echo "Missing/empty raw FASTQ for $sample"
        return 1
    fi


    # Remove incomplete output files before rerunning the sample.
    rm -f "$OUT1" "$OUT2" "$HTML" "$JSON"


    # Run fastp with a time limit and save its exit status.
    set +e

    timeout "$TIMEOUT_FASTP_EXPLICIT" fastp \
        -w "$THREADS" \
        -i "$R1" \
        -I "$R2" \
        -o "$OUT1" \
        -O "$OUT2" \
        --adapter_sequence "$ADAPTER_R1" \
        --adapter_sequence_r2 "$ADAPTER_R2" \
        -h "$HTML" \
        -j "$JSON" \
        > "$LOG" 2>&1

    status=$?

    set -e


    # Move problematic samples to the next retry step.
    if [[ "$status" -ne 0 || ! -s "$OUT1" || ! -s "$OUT2" ]]; then

        echo "FastP explicit failed/timed out for $sample"

        rm -f "$OUT1" "$OUT2"

        mkdir -p "$fail_dir"
        mv "$R1" "$R2" "$fail_dir"/

        return 1
    fi


    # Archive the raw FASTQ files and record the successful run.
    mv "$R1" "$R2" 00_raw_fastq/

    echo -e "${sample}\tfastp_explicit\t${OUT1}\t${OUT2}" >> "$DONE"

    echo "Finished FastP explicit: $sample"

    return 0
}



# Use cutadapt as a fallback for samples that still fail with fastp.
run_cutadapt () {

    sample="$1"
    input_dir="$2"
    fail_dir="$3"

    # Set input, output and log files for this sample.
    R1="${input_dir}/${sample}_R1.fastq.gz"
    R2="${input_dir}/${sample}_R2.fastq.gz"

    OUT1="01_fastp/trimmed_fastq/${sample}_trimmed_R1.fastq.gz"
    OUT2="01_fastp/trimmed_fastq/${sample}_trimmed_R2.fastq.gz"

    LOG="01_fastp/logs/CutadaptLog_${sample}.txt"

    echo "--------------------------------------"
    echo "Cutadapt fallback: $sample"
    echo "Input dir: $input_dir"
    echo "Started: $(date)"


    # Skip samples that were already processed.
    if [[ -s "$OUT1" && -s "$OUT2" ]]; then
        echo "Skipping $sample: trimmed files already exist."
        return 0
    fi


    # Stop this attempt if one of the paired FASTQ files is missing or empty.
    if [[ ! -s "$R1" || ! -s "$R2" ]]; then
        echo "Missing/empty raw FASTQ for $sample"
        return 1
    fi


    # Remove incomplete output files before rerunning the sample.
    rm -f "$OUT1" "$OUT2"


    # Run cutadapt with a time limit and save its exit status.
    set +e

    timeout "$TIMEOUT_CUTADAPT" cutadapt \
        -j "$THREADS" \
        -a "$ADAPTER_R1" \
        -A "$ADAPTER_R2" \
        -o "$OUT1" \
        -p "$OUT2" \
        "$R1" \
        "$R2" \
        > "$LOG" 2>&1

    status=$?

    set -e


    # Move files that still fail to the final failed folder.
    if [[ "$status" -ne 0 || ! -s "$OUT1" || ! -s "$OUT2" ]]; then

        echo "Cutadapt failed/timed out for $sample"

        rm -f "$OUT1" "$OUT2"

        mkdir -p "$fail_dir"
        mv "$R1" "$R2" "$fail_dir"/

        return 1
    fi


    # Archive the raw FASTQ files and record the successful run.
    mv "$R1" "$R2" 00_raw_fastq/

    echo -e "${sample}\tcutadapt\t${OUT1}\t${OUT2}" >> "$DONE"

    echo "Finished Cutadapt fallback: $sample"

    return 0
}



# Run one trimming method on all samples currently in the input folder.
run_round () {

    input_dir="$1"
    method="$2"
    fail_dir="$3"

    echo "======================================"
    echo "Starting round: $method"
    echo "Input dir: $input_dir"
    echo "Fail dir: $fail_dir"
    echo "======================================"


    # Find all R1 files and use their names as sample IDs.
    mapfile -t samples < <(
        find "$input_dir" -maxdepth 1 -name "*_R1.fastq.gz" -printf "%f\n" \
        | sed 's/_R1.fastq.gz//' \
        | sort
    )


    # Skip this round if no samples were found.
    if [[ "${#samples[@]}" -eq 0 ]]; then
        echo "No samples found in $input_dir."
        return 0
    fi


    # Process each sample with the selected trimming method.
    for sample in "${samples[@]}"
    do
        case "$method" in

            fastp_normal)
                run_fastp_normal "$sample" "$input_dir" "$fail_dir" || true
                ;;

            fastp_explicit)
                run_fastp_explicit "$sample" "$input_dir" "$fail_dir" || true
                ;;

            cutadapt)
                run_cutadapt "$sample" "$input_dir" "$fail_dir" || true
                ;;

        esac
    done
}



# Remove old test output files.
rm -f test_*_trimmed_R1.fastq.gz test_*_trimmed_R2.fastq.gz


# First attempt with the standard fastp settings.
run_round "." "fastp_normal" "00_retry_fastp_explicit"


# Retry previously problematic or failed files with the standard fastp settings.
run_round "00_problem_fastq" "fastp_normal" "00_retry_fastp_explicit"
run_round "00_final_failed" "fastp_normal" "00_retry_fastp_explicit"


# Retry remaining problematic files with fixed adapter sequences.
run_round "00_retry_fastp_explicit" "fastp_explicit" "00_retry_cutadapt"


# Use cutadapt for files that still cannot be processed with fastp.
run_round "00_retry_cutadapt" "cutadapt" "00_final_failed"


# Record samples that still have no trimmed paired-end output.
for R1 in 00_final_failed/*_R1.fastq.gz
do
    sample=$(basename "$R1" _R1.fastq.gz)

    if [[ ! -s "01_fastp/trimmed_fastq/${sample}_trimmed_R1.fastq.gz" || \
          ! -s "01_fastp/trimmed_fastq/${sample}_trimmed_R2.fastq.gz" ]]; then

        echo -e "${sample}\tfailed_after_fastp_and_cutadapt" >> "$FAILED"
    fi
done


# Print a short summary of the run.
echo "======================================"
echo "Run finished: $(date)"

echo "Completed samples in this run:"
tail -n +2 "$DONE" | wc -l

echo "Final failed samples:"
tail -n +2 "$FAILED" | wc -l

echo "Completed list: $DONE"
echo "Final failed list: $FAILED"
echo "Main log: $MAIN_LOG"

echo "Total trimmed R1 files now:"
ls 01_fastp/trimmed_fastq/*_trimmed_R1.fastq.gz 2>/dev/null | wc -l

echo "Total trimmed R2 files now:"
ls 01_fastp/trimmed_fastq/*_trimmed_R2.fastq.gz 2>/dev/null | wc -l

echo "======================================"