#!/bin/bash

# This script creates normalized heatmaps and average signal plots for the ChIP-seq data.
# It uses the spike-in normalized ChIP/Input BigWig files to compare H3K9me3 and H4K20me3
# signal between WT and the different KO cell lines.
#
# Peaks from all replicates of the same antibody are first combined into union regions.
# Overlapping peaks are merged, and only regions detected in at least two ChIP replicates
# are kept. This avoids basing the heatmap on peaks that appear in only one replicate.
# deepTools computeMatrix is then used to extract the signal from -5 kb to +5 kb around
# the center of each union region.
#
# Si11 and Si12 are averaged for each cell line, and the regions are sorted by the WT/TT2
# signal. This keeps the same region order for all genotypes and makes differences between
# WT and KO cell lines easier to compare.


# Make sure the correct Conda environment is active.
if [ "${CONDA_DEFAULT_ENV:-}" != "NGS-py37" ]; then
    echo "Error: activate NGS-py37 first." >&2
    exit 1
fi


# Stop the script if a command fails, a variable is missing or a pipeline fails.
set -euo pipefail


# Number of processors used by deepTools computeMatrix.
# Use the number assigned by SLURM if available, otherwise use 32.
THREADS="${SLURM_CPUS_PER_TASK:-32}"
export THREADS

echo "Using ${THREADS} processors for computeMatrix."


# Create the folder used for the peak-signal analysis.
mkdir -p 05_analysis/02_peak_signal/final_overview


# Run the main heatmap analysis in Python.
python3 << 'PY'

import glob          # Find peak files matching the expected file names.
import gzip          # Read the compressed computeMatrix output.
import os            # Work with paths, folders and environment variables.
import shutil        # Find the computeMatrix executable.
import subprocess    # Run deepTools computeMatrix from Python.
import sys           # Stop the script with an error message if required.

import matplotlib.pyplot as plt  # Create the heatmaps and average signal plots.
import numpy as np               # Store, reshape, sort and average signal matrices.


# Analyze 5 kb upstream and downstream of each union-region center
# using 500 bp windows, matching the resolution of the input BigWig files.
BIN = 500
UPSTREAM = 5000
DOWNSTREAM = 5000
NBINS = (UPSTREAM + DOWNSTREAM) // BIN


# Keep regions only if a MACS2 peak was found in at least two ChIP replicates.
MIN_SAMPLE_COUNT = 2


# Use the processor number passed from Bash and make sure that at least one is used.
THREADS = max(1, int(os.environ.get("THREADS", "32")))


# Input and output folders.
BW_DIR = "04_deeptools/bigwig_chip_input_ratio"
PEAK_DIR = "03_macs2/peaks"
OUT_DIR = "05_analysis/02_peak_signal/normalized_heatmaps"

OUTPDF = f"{OUT_DIR}/heatmap_spikein_ratio_union_min2_Si11Si12_merged.pdf"
OUTPNG = f"{OUT_DIR}/heatmap_spikein_ratio_union_min2_Si11Si12_merged.png"

os.makedirs(OUT_DIR, exist_ok=True)


# Antibodies included in the analysis.
MARKS = ["39162", "D84D2"]


# Labels used for the two histone modifications in the final figure.
MARK_LABELS = {
    "39162": "39162 / H3K9me3",
    "D84D2": "D84D2 / H4K20me3",
}


# Order in which the cell lines are shown in the heatmaps and signal plots.
CELL_ORDER = [
    "TT2", "A21", "B23", "G2", "BA3",
    "GB18", "GA6", "A6OHT", "3S1GOHT", "S64-2",
]


# Genotype information shown next to the heatmaps.
KO_INFO = {
    "TT2": "Wildtype",
    "A21": "HP1α KO",
    "B23": "HP1β KO",
    "G2": "HP1γ KO",
    "BA3": "HP1βα DKO",
    "GB18": "HP1γβ DKO",
    "GA6": "HP1γα DKO",
    "A6OHT": "HP1αβγ TKO",
    "3S1GOHT": "H3K9MTs QKO",
    "S64-2": "H4K20MTs DKO",
}


# Colors used for each cell line in the average signal plots and genotype table.
CELL_COLORS = {
    "TT2": "#BFE9E4",
    "A21": "#A8E6CF",
    "B23": "#55B7B1",
    "G2": "#B79AD6",
    "BA3": "#7EA6D9",
    "GB18": "#D65F8C",
    "GA6": "#E8899C",
    "A6OHT": "#D65F8C",
    "3S1GOHT": "#B79AD6",
    "S64-2": "#6C63A8",
}


# Read all MACS2 peak regions belonging to one antibody.
def read_all_peaks_for_mark(mark):

    peak_files = sorted(
        glob.glob(os.path.join(PEAK_DIR, f"*_{mark}_peaks.narrowPeak"))
    )


    # Stop if no peak files were found for this antibody.
    if not peak_files:
        sys.exit(f"Error: no peak files found for {mark}")


    # Store the position, signalValue and replicate name of every peak.
    peaks = []


    # Read the peaks from all replicates belonging to this antibody.
    for peak_file in peak_files:

        sample_name = os.path.basename(peak_file).replace(
            "_peaks.narrowPeak", ""
        )

        with open(peak_file) as handle:

            for line in handle:

                # Ignore empty lines and comment lines in the peak files.
                if not line.strip() or line.startswith("#"):
                    continue


                # Read chromosome, start/end position and MACS2 signalValue.
                fields = line.rstrip().split("\t")

                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                score = float(fields[6]) if len(fields) >= 7 else 0.0

                peaks.append((chrom, start, end, score, sample_name))


    return peaks


# Merge overlapping MACS2 peaks from the different replicates into union regions.
# A union region spans the complete genomic interval covered by overlapping peaks.
def merge_peaks(peaks):

    # Sort the peaks by chromosome and position so overlapping regions are next to each other.
    peaks_sorted = sorted(peaks, key=lambda x: (x[0], x[1], x[2]))

    merged = []


    # Keep track of the union region that is currently being built.
    cur_chr = cur_start = cur_end = None
    cur_score = 0.0
    cur_samples = set()


    # Go through all peaks and combine those that overlap.
    for chrom, start, end, score, sample in peaks_sorted:

        # Start the first union region.
        if cur_chr is None:

            cur_chr, cur_start, cur_end = chrom, start, end
            cur_score = score
            cur_samples = {sample}


        # Extend the current union region when the next peak overlaps it.
        elif chrom == cur_chr and start <= cur_end:

            cur_end = max(cur_end, end)
            cur_score = max(cur_score, score)
            cur_samples.add(sample)


        # The next peak no longer overlaps, so save the finished union region.
        else:

            center = (cur_start + cur_end) // 2

            # Store the genomic range, center, highest signalValue and number
            # of different replicates that contributed a peak to this region.
            merged.append(
                (
                    cur_chr,
                    cur_start,
                    cur_end,
                    center,
                    cur_score,
                    len(cur_samples),
                )
            )


            # Start a new union region with the current non-overlapping peak.
            cur_chr, cur_start, cur_end = chrom, start, end
            cur_score = score
            cur_samples = {sample}


    # Save the last union region after all peaks have been checked.
    if cur_chr is not None:

        center = (cur_start + cur_end) // 2

        merged.append(
            (
                cur_chr,
                cur_start,
                cur_end,
                center,
                cur_score,
                len(cur_samples),
            )
        )


    return merged


# Create union regions by merging overlapping peaks from all replicates of one antibody,
# then keep only regions supported by at least two replicates.
def create_union_regions(mark):

    peaks = read_all_peaks_for_mark(mark)
    merged = merge_peaks(peaks)


    # Remove union regions that contain a peak from only one replicate.
    filtered = [
        region
        for region in merged
        if region[5] >= MIN_SAMPLE_COUNT
    ]


    # BED file used as the region input for deepTools computeMatrix.
    out_bed = os.path.join(
        OUT_DIR,
        f"union_{mark}_min{MIN_SAMPLE_COUNT}.bed"
    )


    # Write the filtered union regions to BED format.
    # The last column records how many replicates contributed a peak to each region.
    with open(out_bed, "w") as out:

        for chrom, start, end, center, score, sample_count in filtered:

            out.write(
                f"{chrom}\t{start}\t{end}\t"
                f"{mark}_union_min{MIN_SAMPLE_COUNT}\t"
                f"{score:.4f}\t.\t{sample_count}\n"
            )


    print(
        f"{mark}: raw peaks = {len(peaks)}, merged regions = {len(merged)}, "
        f"kept regions with >= {MIN_SAMPLE_COUNT} peaksets = {len(filtered)}"
    )


    # Return the chromosome and center used as reference points for the signal analysis.
    return [
        (chrom, center)
        for chrom, start, end, center, score, sample_count in filtered
    ]


# Build the expected BigWig path for one cell line, replicate and antibody.
def bw_path(cell, si, mark):

    sample = f"{cell}-{si}_{mark}"

    return os.path.join(
        BW_DIR,
        f"{sample}_vs_Inp_spikein_ratio_bin{BIN}.bw"
    )


# Extract the normalized BigWig signal around the union regions with computeMatrix.
def compute_matrix_for_mark(mark, regions):

    # Find the deepTools computeMatrix program in the active environment.
    compute_matrix_exe = shutil.which("computeMatrix")


    # Stop if computeMatrix is not available.
    if compute_matrix_exe is None:

        sys.exit(
            "Error: computeMatrix was not found in PATH. "
            "Please activate the environment that contains deepTools."
        )


    # BED file containing the union regions for this antibody.
    region_bed = os.path.join(
        OUT_DIR,
        f"union_{mark}_min{MIN_SAMPLE_COUNT}.bed"
    )


    bigwigs = []
    sample_labels = []


    # Collect Si11 and Si12 for every cell line in a fixed order.
    # The same order is later used to average the two replicates correctly.
    for cell in CELL_ORDER:

        for si in ("Si11", "Si12"):

            path = bw_path(cell, si, mark)


            # Stop if one of the BigWigs needed for the comparison is missing.
            if not os.path.exists(path):
                sys.exit(f"Error: missing BigWig: {path}")


            bigwigs.append(path)
            sample_labels.append(f"{cell}-{si}")


    # Output files created by computeMatrix.
    matrix_gz = os.path.join(
        OUT_DIR,
        f"computeMatrix_{mark}_union_min{MIN_SAMPLE_COUNT}.gz"
    )

    matrix_tab = os.path.join(
        OUT_DIR,
        f"computeMatrix_{mark}_union_min{MIN_SAMPLE_COUNT}.tab"
    )


    # Extract the signal from all BigWigs in 500 bp windows around
    # the center of every union region.
    cmd = [
        compute_matrix_exe,
        "reference-point",
        "--referencePoint", "center",
        "-R", region_bed,
        "-S", *bigwigs,
        "-b", str(UPSTREAM),
        "-a", str(DOWNSTREAM),
        "--binSize", str(BIN),
        "--missingDataAsZero",
        "--numberOfProcessors", str(THREADS),
        "--samplesLabel", *sample_labels,
        "-o", matrix_gz,
        "--outFileNameMatrix", matrix_tab,
    ]


    print("", flush=True)
    print("=" * 72, flush=True)

    print(
        f"{mark}: computeMatrix starting with {THREADS} processors "
        f"for {len(bigwigs)} BigWigs and {len(regions):,} regions",
        flush=True
    )

    print("=" * 72, flush=True)


    # Run deepTools computeMatrix.
    subprocess.run(cmd, check=True)


    print(
        f"{mark}: computeMatrix finished.",
        flush=True
    )

    print(
        f"{mark}: loading matrix and averaging Si11/Si12...",
        flush=True
    )


    # Each BigWig contributes one value for every 500 bp window.
    expected_signal_columns = len(bigwigs) * NBINS

    rows = []


    # Read the signal values from the compressed computeMatrix output.
    with gzip.open(matrix_gz, "rt") as handle:

        for line in handle:

            # Ignore empty lines and computeMatrix header lines.
            if not line.strip() or line.startswith("@"):
                continue


            fields = line.rstrip("\n").split("\t")


            # Stop if a row does not contain all expected signal values.
            if len(fields) < expected_signal_columns:

                sys.exit(
                    f"Error while parsing {matrix_gz}: row has too few columns."
                )


            # Keep only the columns containing BigWig signal values.
            signal_fields = fields[-expected_signal_columns:]


            # Convert the signal values to numbers and use zero where no value is available.
            values = np.array(
                [
                    float(x)
                    if x not in ("nan", "NaN", "NA", "")
                    else 0.0
                    for x in signal_fields
                ],
                dtype=float
            )

            values[~np.isfinite(values)] = 0.0

            rows.append(values)


    # Stop if no region signal was read from the matrix.
    if not rows:
        sys.exit(f"Error: no matrix rows were parsed for {mark}")


    # Reshape the computeMatrix output so the signal can be handled separately
    # for each union region, replicate and 500 bp window.
    raw_matrix = np.vstack(rows).reshape(
        len(rows),
        len(bigwigs),
        NBINS
    )


    merged_cells = []


    # Average Si11 and Si12 so each cell line is represented by one signal matrix.
    for i, cell in enumerate(CELL_ORDER):

        si11_idx = i * 2
        si12_idx = i * 2 + 1

        merged = (
            raw_matrix[:, si11_idx, :] +
            raw_matrix[:, si12_idx, :]
        ) / 2.0

        merged_cells.append(merged)


    # Place the cell-line matrices next to each other for the final heatmap.
    matrix = np.concatenate(merged_cells, axis=1)


    print(
        f"{mark}: replicate averaging finished; matrix shape = {matrix.shape}",
        flush=True
    )


    return matrix


# Build the union-region set separately for H3K9me3 and H4K20me3.
regions_by_mark = {}

for mark in MARKS:

    print(
        f"{mark}: building merged peak regions...",
        flush=True
    )

    regions_by_mark[mark] = create_union_regions(mark)


# Store the final signal matrix for each histone modification.
panel_data = {}


# Build the signal matrix separately for H3K9me3 and H4K20me3,
# using the union regions created for each antibody.
for mark in MARKS:

    matrix = compute_matrix_for_mark(
        mark,
        regions_by_mark[mark]
    )


    # Use the average WT/TT2 signal across the 10 kb region to sort the heatmap rows.
    wt_matrix = matrix[:, 0:NBINS]

    wt_mean_signal = np.nanmean(
        wt_matrix,
        axis=1
    )

    order = np.argsort(-wt_mean_signal)


    # Apply the WT-based row order to all genotypes so the same regions can be compared.
    panel_data[mark] = matrix[order]


    print(
        f"{mark}: WT-based sorting finished.",
        flush=True
    )


print(
    "Signal extraction complete. Creating final figure...",
    flush=True
)


# Create the main figure.
fig = plt.figure(
    figsize=(14.5, 12.5)
)


# Layout for average signal plots, heatmaps, color bars and genotype information.
gs = fig.add_gridspec(
    nrows=4,
    ncols=3,
    width_ratios=[1.0, 0.045, 0.55],
    height_ratios=[0.28, 1.0, 0.28, 1.0],
    wspace=0.16,
    hspace=0.34,
)


# Axes for the average signal plots shown above each heatmap.
profile_axes = {
    "39162": fig.add_subplot(gs[0, 0]),
    "D84D2": fig.add_subplot(gs[2, 0]),
}


# Axes for the two heatmaps.
heat_axes = {
    "39162": fig.add_subplot(gs[1, 0]),
    "D84D2": fig.add_subplot(gs[3, 0]),
}


# Separate signal scale for each histone modification.
cbar_axes = {
    "39162": fig.add_subplot(gs[1, 1]),
    "D84D2": fig.add_subplot(gs[3, 1]),
}


# Area used for the cell-line and genotype information.
legend_ax = fig.add_subplot(gs[:, 2])
legend_ax.axis("off")


# Create the signal plot and heatmap for each histone modification.
for mark in MARKS:

    matrix = panel_data[mark]

    ax_prof = profile_axes[mark]
    ax_heat = heat_axes[mark]

    n_samples = len(CELL_ORDER)


    # Calculate the average signal across all union regions for each cell line.
    for i, cell in enumerate(CELL_ORDER):

        start = i * NBINS
        end = (i + 1) * NBINS

        mean_profile = np.nanmean(
            matrix[:, start:end],
            axis=0
        )


        # Plot how the average signal changes from -5 kb to +5 kb around the region center.
        ax_prof.plot(
            mean_profile,
            color=CELL_COLORS[cell],
            linewidth=1.6,
            alpha=0.95,
        )


    # Add the antibody and corresponding histone modification above the signal plot.
    ax_prof.set_title(
        MARK_LABELS[mark],
        fontsize=13,
        fontweight="bold",
        pad=7
    )


    # Show the genomic position relative to the union-region center.
    ax_prof.set_xlim(0, NBINS - 1)

    ax_prof.set_xticks(
        [0, NBINS // 2, NBINS - 1]
    )

    ax_prof.set_xticklabels(
        ["-5 kb", "center", "+5 kb"],
        fontsize=8
    )

    ax_prof.set_ylabel(
        "Normalized ChIP signal",
        fontsize=8.5,
        labelpad=8
    )

    ax_prof.tick_params(
        axis="y",
        length=0,
        labelleft=False
    )

    ax_prof.tick_params(
        axis="x",
        length=0,
        pad=2
    )


    # Remove the surrounding frame from the signal plot.
    for spine in ax_prof.spines.values():
        spine.set_visible(False)


    # Use the central 98% of signal values for the color scale so a few
    # extreme values do not determine the appearance of the whole heatmap.
    vmin = np.nanpercentile(
        matrix.flatten(),
        1
    )

    vmax = np.nanpercentile(
        matrix.flatten(),
        99
    )


    # Plot the signal for all union regions and cell lines as a heatmap.
    im = ax_heat.imshow(
        matrix,
        aspect="auto",
        interpolation="nearest",
        cmap="inferno",
        vmin=vmin,
        vmax=vmax,
    )


    # Add white lines between the different cell-line blocks.
    for i in range(1, n_samples):

        ax_heat.axvline(
            i * NBINS - 0.5,
            color="white",
            linewidth=0.6,
            alpha=0.75,
        )


    # Calculate the middle of each cell-line block for the x-axis labels.
    centers = [
        (i * NBINS) + (NBINS / 2) - 0.5
        for i in range(n_samples)
    ]


    # Place each cell-line name below its corresponding heatmap block.
    ax_heat.set_xticks(centers)

    ax_heat.set_xticklabels(
        CELL_ORDER,
        rotation=45,
        ha="right",
        fontsize=8
    )

    ax_heat.tick_params(
        axis="y",
        length=0,
        labelleft=False
    )

    ax_heat.tick_params(
        axis="x",
        length=0,
        pad=1
    )


    # Remove the frame around the heatmap.
    for spine in ax_heat.spines.values():
        spine.set_visible(False)


    ax_heat.set_ylabel("")


    # Add the signal scale next to the corresponding heatmap.
    cbar = fig.colorbar(
        im,
        cax=cbar_axes[mark]
    )

    cbar.set_label(
        "Normalized ChIP signal",
        fontsize=9.2,
        labelpad=10
    )

    cbar.ax.tick_params(
        labelsize=8,
        length=0
    )

    cbar.outline.set_visible(False)


# Add the heading for the cell-line and genotype information.
legend_ax.text(
    0.56,
    0.94,
    "Cell line / Genotype",
    ha="center",
    va="center",
    fontsize=11.5,
    fontweight="bold",
)


# Build the table containing the cell-line and KO information.
table_data = [
    ["Cell line", "KO type"]
]

for cell in CELL_ORDER:
    table_data.append(
        [cell, KO_INFO[cell]]
    )


# Use the same cell-line colors in the table as in the signal plots.
cell_colours = [
    ["#F5F5F5", "#F5F5F5"]
]

for cell in CELL_ORDER:
    cell_colours.append(
        [CELL_COLORS[cell], "white"]
    )


# Add the genotype table to the right side of the figure.
table = legend_ax.table(
    cellText=table_data,
    cellColours=cell_colours,
    cellLoc="center",
    colWidths=[0.42, 0.58],
    bbox=[0.08, 0.34, 0.90, 0.54],
)


# Set the table font size manually.
table.auto_set_font_size(False)
table.set_fontsize(8.5)


# Format the header, cell-line and genotype columns separately.
for (row, col), cell in table.get_celld().items():

    cell.set_edgecolor("none")
    cell.set_linewidth(0)


    if row == 0:

        cell.set_text_props(
            weight="bold",
            color="black"
        )

        cell.set_facecolor("#F5F5F5")


    elif col == 0:

        cell.set_text_props(
            weight="bold",
            color="white"
        )

        cell.set_facecolor(
            cell_colours[row][0]
        )


    else:

        cell.set_text_props(
            color="black"
        )

        cell.set_facecolor("white")


# Add the main title above the complete figure.
fig.suptitle(
    "Normalized ChIP signal heatmap",
    fontsize=16,
    fontweight="bold",
    y=0.985,
)


# Summarize how replicates, union regions and row sorting were handled.
fig.text(
    0.075,
    0.025,
    f"Si11 and Si12 were averaged per cell line. All samples are plotted over the same merged peak regions per antibody. Regions kept if detected in at least {MIN_SAMPLE_COUNT} peaksets. Rows are sorted by average WT/TT2 signal.",
    fontsize=9,
    ha="left",
)


# Save the final figure as both PDF and PNG.
plt.savefig(
    OUTPDF,
    format="pdf",
    dpi=300,
    bbox_inches="tight"
)

plt.savefig(
    OUTPNG,
    format="png",
    dpi=300,
    bbox_inches="tight"
)


# Close the figure after saving.
plt.close()


print(f"Finished: {OUTPDF}")
print(f"Finished: {OUTPNG}")

PY