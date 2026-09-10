#!/bin/bash

# This script identifies enriched ChIP-seq regions with MACS2.
# Each ChIP sample is compared with the matching Input control from the same cell line.
# The Input control contains DNA from the same sample preparation without antibody
# enrichment and is used to estimate background signal.
# MACS2 uses this comparison to distinguish enriched ChIP regions from background.
# Several samples can be processed in parallel to reduce the total runtime.
# The resulting peak files contain genomic regions with significant ChIP enrichment
# and are used for downstream comparison and visualization.


# Stop the script on errors and unset variables.
set -Eeuo pipefail


# Make sure the correct Conda environment is active.
if [[ "${CONDA_DEFAULT_ENV:-}" != "NGS-py27" ]]; then
    echo "ERROR: Please activate NGS-py27 before running this script."
    exit 1
fi


# Maximum number of MACS2 samples processed at the same time.
MAX_JOBS=6


# Check that MACS2 is available.
if ! command -v macs2 >/dev/null 2>&1; then
    echo "Error: macs2 not found in NGS-py27."
    exit 1
fi


echo "Running MACS2 in environment: $CONDA_DEFAULT_ENV"
echo "Parallel MACS2 jobs: $MAX_JOBS"

macs2 --version


# Create output and log folders.
mkdir -p 03_macs2
mkdir -p 03_macs2/logs


# Find all Input BAM files used as background controls for peak calling.
control_bams=$(ls 02_alignment/bam_mm10/*Inp*_mm10_sorted_markdup.bam 2>/dev/null || true)


# Stop if no Input/control BAM files were found.
if [ -z "$control_bams" ]; then
    echo "Error: No Input/control mm10 markdup BAM found in 02_alignment/bam_mm10/"
    exit 1
fi


# Process each Input control separately.
for control_bam in $control_bams
do

    filename=$(basename "$control_bam")

    # Get the cell line name from the Input filename.
    cell=${filename%%_Inp*}


    # Find all ChIP BAM files belonging to the same cell line.
    treatment_bams=$(ls 02_alignment/bam_mm10/"${cell}"*_mm10_sorted_markdup.bam 2>/dev/null | grep -v 'Inp' || true)


    # Skip the cell line if no matching ChIP samples were found.
    if [ -z "$treatment_bams" ]; then
        echo "Warning: No treatment mm10 BAM files found for $cell. Skipping."
        continue
    fi


    # Run MACS2 for each ChIP sample belonging to this Input control.
    for treatment_bam in $treatment_bams
    do

        treatment_filename=$(basename "$treatment_bam")

        # Get the sample name from the BAM filename.
        sample=$(echo "$treatment_filename" | sed -e 's/_mm10_sorted_markdup.bam//g')


        # Define the expected peak file and log file for this sample.
        expected_peak="03_macs2/${sample}_peaks.narrowPeak"
        log_file="03_macs2/logs/MACS2_log_${sample}.txt"


        # Skip samples if the MACS2 peak file already exists.
        if [[ -s "$expected_peak" ]]; then
            echo "Skipping $sample: MACS2 peak file already exists."
            continue
        fi


        echo "Starting MACS2 for $sample"


        # Run MACS2 peak calling using the matching Input as background control.
        (
            macs2 callpeak \
                -t "$treatment_bam" \
                -c "$control_bam" \
                -f BAMPE \
                -g mm \
                -n "$sample" \
                --outdir 03_macs2 \
                -q 0.01 \
                > "$log_file" 2>&1


            # Check that MACS2 created the expected peak file.
            if [[ ! -s "$expected_peak" ]]; then
                echo "Error: MACS2 did not create expected peak file for $sample" >> "$log_file"
                exit 1
            fi


            echo "MACS2 finished for $sample" >> "$log_file"

        ) &


        # Wait if the maximum number of parallel MACS2 jobs is already running.
        while [ "$(jobs -r | wc -l)" -ge "$MAX_JOBS" ]
        do
            sleep 30
        done

    done

done


# Wait until all remaining background MACS2 jobs are finished.
wait


echo "All MACS2 jobs finished."