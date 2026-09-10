#!/bin/bash

# This script creates spike-in normalized BigWig files from the aligned ChIP-seq data.
# Drosophila (dm6) reads are used as the spike-in reference to correct for differences
# between samples. This makes the ChIP signal more comparable across samples even if
# sequencing depth or global signal levels differ.
# First, the number of dm6 reads is counted for each sample and a scaling factor is
# calculated relative to the mean dm6 read count. The scaling factor is then applied
# to the corresponding mouse (mm10) BAM file with bamCoverage.
# Finally, each ChIP BigWig is divided by its matching Input BigWig to generate
# ChIP/Input ratio tracks for downstream visualization and analysis.


# Stop the script on errors and unset variables.
set -Eeuo pipefail

# Ignore file patterns if no matching files are found.
shopt -s nullglob


# Make sure the correct Conda environment is active.
if [[ "${CONDA_DEFAULT_ENV:-}" != "NGS-py37" ]]; then
    echo "ERROR: Please activate NGS-py37 before running this script."
    exit 1
fi


# Number of threads used by deepTools.
THREADS=60

# Size of the genomic windows used to calculate the BigWig signal.
BIN=500


# Create folders for logs, BigWig files and scaling factors.
mkdir -p 04_deeptools/logs
mkdir -p 04_deeptools/bigwig_spikein_scaled
mkdir -p 04_deeptools/bigwig_chip_input_ratio
mkdir -p 04_deeptools/scaling_factors


# Check that all required programs are available.
for cmd in bamCoverage bigwigCompare samtools bc; do
    command -v "$cmd" >/dev/null 2>&1 || {
        echo "ERROR: $cmd not found"
        exit 1
    }
done


# Save all terminal output in a run-specific log file.
LOG_FILE="04_deeptools/logs/deeptools_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1


echo "Running deepTools in: $CONDA_DEFAULT_ENV"
echo "Threads: $THREADS"
echo "Bin size: $BIN"


# File used to store the dm6 read counts and calculated scaling factors.
SCALING_FILE="04_deeptools/scaling_factors/dm6_scaling_factors.tsv"
echo -e "sample\tdm6_reads\tscale_factor" > "$SCALING_FILE"


# Count mapped dm6 reads for each sample.
DM6_COUNTS="04_deeptools/scaling_factors/dm6_read_counts.tsv"
echo -e "sample\tdm6_reads" > "$DM6_COUNTS"

for dm6_bam in 02_alignment/bam_dm6/*_dm6_sorted_markdup.bam; do

    sample=$(basename "$dm6_bam" _dm6_sorted_markdup.bam)

    # Count mapped Drosophila reads after duplicate removal.
    reads=$(samtools view -c -F 4 "$dm6_bam")

    echo -e "${sample}\t${reads}" >> "$DM6_COUNTS"

done


# Calculate the mean dm6 read count across all samples.
# This is used as the reference for the spike-in scaling factors.
MEAN_DM6=$(awk 'NR>1 {sum+=$2; n++} END {if(n>0) printf "%.10f\n", sum/n; else print 0}' "$DM6_COUNTS")

if (( $(echo "$MEAN_DM6 <= 0" | bc -l) )); then
    echo "ERROR: Mean dm6 read count is 0."
    exit 1
fi

echo "Mean dm6 reads: $MEAN_DM6"


# Calculate the spike-in scaling factor for each sample
# and create the scaled mm10 BigWig file.
for mm10_bam in 02_alignment/bam_mm10/*_mm10_sorted_markdup.bam; do

    sample=$(basename "$mm10_bam" _mm10_sorted_markdup.bam)

    # Get the dm6 read count belonging to the same sample.
    dm6_reads=$(awk -v s="$sample" '$1==s {print $2}' "$DM6_COUNTS")


    # Stop if the matching dm6 read count is missing or zero.
    if [[ -z "$dm6_reads" || "$dm6_reads" -eq 0 ]]; then
        echo "ERROR: Missing/zero dm6 reads for $sample"
        exit 1
    fi


    # Samples with fewer dm6 reads receive a larger scaling factor and vice versa.
    scale_factor=$(awk -v mean="$MEAN_DM6" -v reads="$dm6_reads" \
        'BEGIN {printf "%.10f\n", mean/reads}')

    echo -e "${sample}\t${dm6_reads}\t${scale_factor}" >> "$SCALING_FILE"


    out_bw="04_deeptools/bigwig_spikein_scaled/${sample}_spikein_scaled_bin${BIN}.bw"


    # Skip the sample if the scaled BigWig already exists.
    if [[ -s "$out_bw" ]]; then
        echo "Skipping existing scaled BigWig: $sample"
        continue
    fi


    echo "Creating spike-in scaled BigWig for $sample"
    echo "  dm6 reads: $dm6_reads"
    echo "  scale factor: $scale_factor"


    # Convert the mm10 BAM file to BigWig and apply the dm6 scaling factor.
    bamCoverage \
        -b "$mm10_bam" \
        -o "$out_bw" \
        -bs "$BIN" \
        -p "$THREADS" \
        --scaleFactor "$scale_factor"

done


# Find Input BigWigs and match them to ChIP samples from the same cell line.
for input_bw in 04_deeptools/bigwig_spikein_scaled/*Inp*_spikein_scaled_bin${BIN}.bw; do

    input_file=$(basename "$input_bw")
    cell=${input_file%%_Inp*}


    # Find the ChIP BigWigs belonging to the same cell line.
    for chip_bw in 04_deeptools/bigwig_spikein_scaled/${cell}*_spikein_scaled_bin${BIN}.bw; do

        # Do not compare Input samples against themselves.
        [[ "$chip_bw" == *Inp* ]] && continue


        # Remove the path and suffix from the ChIP BigWig name and define the ratio output file.
        chip_file=$(basename "$chip_bw")
        sample=${chip_file%_spikein_scaled_bin${BIN}.bw}

        out_ratio="04_deeptools/bigwig_chip_input_ratio/${sample}_vs_Inp_spikein_ratio_bin${BIN}.bw"


        # Skip the comparison if the ratio BigWig already exists.
        if [[ -s "$out_ratio" ]]; then
            echo "Skipping existing ChIP/Input ratio: $sample"
            continue
        fi


        echo "Creating ChIP/Input ratio for $sample"


        # Divide the spike-in scaled ChIP signal by the matching Input signal.
        bigwigCompare \
            -b1 "$chip_bw" \
            -b2 "$input_bw" \
            -o "$out_ratio" \
            --operation ratio \
            -bs "$BIN" \
            -p "$THREADS"

    done

done


echo "deepTools processing completed."
echo "Scaling factors: $SCALING_FILE"
echo "Log: $LOG_FILE"