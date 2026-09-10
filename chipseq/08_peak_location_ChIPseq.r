# This script summarizes the genomic locations of the annotated ChIP-seq peaks
# and creates plots showing both absolute peak numbers and relative percentages.
#
# It uses the annotation tables created by the previous ChIPseeker script and
# groups the detailed ChIPseeker annotations into broader genomic categories
# such as promoter, exon, intron, UTR, downstream and intergenic regions.
#
# Peak numbers are first calculated separately for Si11 and Si12 and then averaged
# for each cell line. This allows the genomic distribution of H3K9me3 and H4K20me3
# peaks to be compared across WT and the different KO genotypes.
#
# Two plot types are created: one shows the average number of peaks in each genomic
# category, while the second shows the same distribution as percentages so that
# relative differences can be compared independently of the total number of peaks.


# Defines the CRAN repository used if an R package has to be installed.
# CRAN is the main online repository for general R packages.
cran_repo <- "https://cloud.r-project.org"


# Lists the R packages used for data processing and plotting.
required_packages <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr",
  "readr"
)


# Checks whether each required package is available and installs missing packages.
for (pkg in required_packages) {

  if (!requireNamespace(pkg, quietly = TRUE)) {

    install.packages(
      pkg,
      repos = cran_repo
    )
  }
}


# dplyr is used to filter, group, count and modify the annotation data.
library(dplyr)

# tidyr is used to add missing rows for sample, replicate and genomic-region groups.
library(tidyr)

# ggplot2 is used to create the absolute and percentage bar plots.
library(ggplot2)

# stringr is used to identify information in filenames and annotation text.
library(stringr)

# readr provides functions for reading and writing tabular data.
library(readr)


# -----------------------
# Paths
# -----------------------

# Folder containing the peak annotation tables created by the previous script.
annotation_dir <- "05_analysis/04_annotation/annotation_tables"

# Folder where the summarized tables and plots are saved.
out_dir <- "05_analysis/04_annotation/location_plots"


# Creates the output folder, including parent folders if they do not exist yet.
dir.create(
  out_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# Defines the output files for the absolute peak numbers and percentages.
out_counts_csv <- file.path(
  out_dir,
  "peak_counts_by_genomic_location.csv"
)

out_percent_csv <- file.path(
  out_dir,
  "peak_percent_by_genomic_location.csv"
)

out_counts_pdf <- file.path(
  out_dir,
  "peak_counts_by_genomic_location.pdf"
)

out_counts_png <- file.path(
  out_dir,
  "peak_counts_by_genomic_location.png"
)

out_percent_pdf <- file.path(
  out_dir,
  "peak_percent_by_genomic_location.pdf"
)

out_percent_png <- file.path(
  out_dir,
  "peak_percent_by_genomic_location.png"
)


# -----------------------
# Settings
# -----------------------

# Defines the order in which the cell lines are shown in the plots.
cell_order <- c(
  "TT2",
  "A21",
  "B23",
  "G2",
  "BA3",
  "GB18",
  "GA6",
  "A6OHT",
  "3S1GOHT",
  "S64-2"
)


# Defines the labels used for the two antibodies/histone modifications.
mark_labels <- c(
  "39162" = "H3K9me3 / 39162",
  "D84D2" = "H4K20me3 / D84D2"
)


# Defines the genomic categories and their order in the plots.
# UTR means untranslated region, which is part of an RNA transcript but is
# located outside the protein-coding sequence.
region_order <- c(
  "Promoter",
  "5UTR",
  "3UTR",
  "Exon",
  "Intron",
  "Downstream",
  "Intergenic",
  "Other"
)


# Defines the colors used for the different genomic-location categories.
region_colors <- c(
  "Promoter" = "#D65F8C",
  "5UTR" = "#E8899C",
  "3UTR" = "#F2B880",
  "Exon" = "#A8E6CF",
  "Intron" = "#7EA6D9",
  "Downstream" = "#55B7B1",
  "Intergenic" = "#B79AD6",
  "Other" = "#BDBDBD"
)


# -----------------------
# Helper functions
# -----------------------

# Converts the detailed ChIPseeker annotation into broader genomic categories.
classify_region <- function(annotation) {

  # Converts the annotation text to lower case so capitalization does not
  # affect the category matching.
  annotation <- tolower(annotation)


  # Checks the annotation text for characteristic terms and assigns each peak
  # to one of the genomic categories used in the plots.
  case_when(
    str_detect(annotation, "promoter") ~ "Promoter",
    str_detect(annotation, "5.?utr|five.?utr") ~ "5UTR",
    str_detect(annotation, "3.?utr|three.?utr") ~ "3UTR",
    str_detect(annotation, "exon") ~ "Exon",
    str_detect(annotation, "intron") ~ "Intron",
    str_detect(annotation, "downstream") ~ "Downstream",
    str_detect(annotation, "distal intergenic|intergenic") ~ "Intergenic",
    TRUE ~ "Other"
  )
}


# Extracts the cell line, replicate and antibody from the annotation filename.
parse_file_info <- function(file) {

  # Removes the folder path and keeps only the filename.
  fname <- basename(file)


  # Separates filenames such as TT2-Si11_39162_peak_annotation.csv
  # into cell line, replicate and antibody information.
  parsed <- str_match(
    fname,
    "^(.+)-(Si[0-9]+)_(39162|D84D2)_peak_annotation\\.csv$"
  )


  # Stops the script if the filename does not follow the expected naming pattern,
  # because the sample information could otherwise be assigned incorrectly.
  if (any(is.na(parsed))) {

    stop(
      paste(
        "Could not parse annotation filename:",
        fname
      )
    )
  }


  # Returns the extracted sample information as a small data frame.
  data.frame(
    file = file,
    filename = fname,
    cell = parsed[, 2],
    si = parsed[, 3],
    mark = parsed[, 4],
    stringsAsFactors = FALSE
  )
}


# Reads one annotation table and prepares the information needed for the plots.
read_annotation_file <- function(file) {

  # Extracts cell line, replicate and antibody information from the filename.
  info <- parse_file_info(file)


  # Reads the complete ChIPseeker annotation table.
  df <- read.csv(
    file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )


  # Finds the column containing the genomic annotation.
  # grep searches the column names for either "annotation" or "annot".
  annotation_col <- grep(
    "annotation|annot",
    colnames(df),
    value = TRUE,
    ignore.case = TRUE
  )[1]


  # Stops the script if no annotation column can be identified.
  if (is.na(annotation_col)) {

    stop(
      paste(
        "No annotation column found in:",
        file
      )
    )
  }


  # Finds the column containing the chromosome information.
  chr_col <- grep(
    "^seqnames$|^chr$|chrom",
    colnames(df),
    value = TRUE,
    ignore.case = TRUE
  )[1]


  # Finds the column containing the peak start position.
  start_col <- grep(
    "^start$",
    colnames(df),
    value = TRUE,
    ignore.case = TRUE
  )[1]


  # Finds the column containing the peak end position.
  end_col <- grep(
    "^end$",
    colnames(df),
    value = TRUE,
    ignore.case = TRUE
  )[1]


  # Creates a unique peak ID from chromosome, start and end coordinates when
  # these columns are available. Otherwise, the row number is used as the ID.
  if (
    !is.na(chr_col) &&
    !is.na(start_col) &&
    !is.na(end_col)
  ) {

    df$peak_id <- paste(
      df[[chr_col]],
      df[[start_col]],
      df[[end_col]],
      sep = "_"
    )

  } else {

    df$peak_id <- seq_len(
      nrow(df)
    )
  }


  # Adds the sample information and converts the detailed ChIPseeker annotation
  # into the broader genomic categories used for the final plots.
  # Duplicate entries for the same peak and category are removed.
  df %>%
    mutate(
      cell = info$cell,
      si = info$si,
      mark = info$mark,
      region = classify_region(
        .data[[annotation_col]]
      )
    ) %>%
    distinct(
      cell,
      si,
      mark,
      peak_id,
      region
    )
}


# -----------------------
# Read files
# -----------------------

# Finds all peak annotation tables created by the previous annotation script.
annotation_files <- list.files(
  annotation_dir,
  pattern = "_peak_annotation\\.csv$",
  full.names = TRUE
)


# Stops the script if no annotation tables were found.
if (length(annotation_files) == 0) {

  stop(
    paste(
      "No peak annotation files found in:",
      annotation_dir
    )
  )
}


# Prints the number of annotation files found.
message(
  "Found annotation files: ",
  length(annotation_files)
)


# Reads every annotation table and combines all samples into one data frame.
all_annotations <- bind_rows(
  lapply(
    annotation_files,
    read_annotation_file
  )
)


# Keeps only the cell lines and antibodies included in this analysis.
all_annotations <- all_annotations %>%
  filter(
    cell %in% cell_order,
    mark %in% names(mark_labels)
  )


# -----------------------
# Count peaks per replicate
# -----------------------

# Counts the number of peaks in each genomic category separately for every
# antibody, cell line and replicate.
rep_counts <- all_annotations %>%
  count(
    mark,
    cell,
    si,
    region,
    name = "peak_count"
  ) %>%

  # Adds rows for genomic categories with no detected peaks and assigns
  # them a count of zero so every sample contains the same categories.
  complete(
    mark,
    cell = cell_order,
    si,
    region = region_order,
    fill = list(
      peak_count = 0
    )
  )


# -----------------------
# Average Si11 and Si12 per cell line
# -----------------------

# Averages the peak counts of Si11 and Si12 so each cell line is represented
# by one value for every antibody and genomic category.
avg_counts <- rep_counts %>%
  group_by(
    mark,
    cell,
    region
  ) %>%
  summarise(
    mean_peak_count = mean(
      peak_count,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%

  # Converts the variables to factors so the plots follow the defined order
  # instead of using alphabetical order.
  mutate(
    cell = factor(
      cell,
      levels = cell_order
    ),
    region = factor(
      region,
      levels = region_order
    ),
    mark_label = factor(
      mark_labels[mark],
      levels = mark_labels[
        c(
          "39162",
          "D84D2"
        )
      ]
    )
  )


# -----------------------
# Percent per cell line and mark
# -----------------------

# Calculates the relative percentage of peaks in each genomic category
# within every cell line and antibody.
percent_counts <- avg_counts %>%
  group_by(
    mark,
    mark_label,
    cell
  ) %>%
  mutate(

    # Calculates the total average peak number across all genomic categories.
    total_peaks = sum(
      mean_peak_count,
      na.rm = TRUE
    ),

    # Calculates the percentage contributed by each genomic category.
    # A value of zero is used if no peaks are present.
    percent_peaks = ifelse(
      total_peaks > 0,
      100 * mean_peak_count / total_peaks,
      0
    )
  ) %>%
  ungroup()


# Saves the absolute and percentage summaries as CSV tables.
write.csv(
  avg_counts,
  out_counts_csv,
  row.names = FALSE
)

write.csv(
  percent_counts,
  out_percent_csv,
  row.names = FALSE
)


# -----------------------
# Shared plot theme
# -----------------------

# Defines the common appearance used for both plots so their layout and
# formatting remain consistent.
base_theme <- theme_minimal(
  base_size = 12
) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 10
    ),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      face = "bold",
      size = 9
    ),
    axis.title.x = element_blank(),
    axis.title.y = element_text(
      face = "bold"
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    strip.background = element_rect(
      fill = "#F1F1F1",
      color = NA
    ),
    legend.position = "right",
    legend.title = element_text(
      face = "bold"
    )
  )


# -----------------------
# Plot 1: absolute peak numbers
# -----------------------

# Creates stacked bar plots showing the average number of peaks in each
# genomic category for every cell line.
p_counts <- ggplot(
  avg_counts,
  aes(
    x = cell,
    y = mean_peak_count,
    fill = region
  )
) +

  # Draws one stacked bar per cell line, with one section for each genomic category.
  geom_col(
    width = 0.68,
    color = "white",
    linewidth = 0.2
  ) +

  # Shows H3K9me3 and H4K20me3 in separate panels because their total peak
  # numbers can differ and therefore use separate y-axis ranges.
  facet_wrap(
    ~ mark_label,
    ncol = 1,
    scales = "free_y"
  ) +

  # Uses the previously defined color for each genomic category.
  scale_fill_manual(
    values = region_colors,
    drop = FALSE
  ) +

  base_theme +

  # Adds the plot title, subtitle, y-axis label and legend title.
  labs(
    title = "Number of H3K9me3 and H4K20me3 peaks by genomic location",
    subtitle = "Si11 and Si12 peak counts were averaged per cell line",
    y = "Number of peaks",
    fill = "Genomic location"
  )


# -----------------------
# Plot 2: percentage distribution
# -----------------------

# Creates stacked percentage plots showing how the peaks of each cell line
# are distributed across the different genomic categories.
p_percent <- ggplot(
  avg_counts,
  aes(
    x = cell,
    y = mean_peak_count,
    fill = region
  )
) +

  # Scales each complete bar to 100%, so the different genomic categories
  # represent their relative contribution instead of absolute peak numbers.
  geom_col(
    width = 0.68,
    color = "white",
    linewidth = 0.2,
    position = "fill"
  ) +

  # Shows H3K9me3 and H4K20me3 in separate panels.
  facet_wrap(
    ~ mark_label,
    ncol = 1
  ) +

  scale_fill_manual(
    values = region_colors,
    drop = FALSE
  ) +

  # Converts the internal 0-1 scale created by position = "fill"
  # into percentage labels from 0 to 100.
  scale_y_continuous(
    breaks = seq(
      0,
      1,
      0.2
    ),
    labels = function(x) x * 100
  ) +

  base_theme +

  # Adds the plot title, subtitle, percentage y-axis label and legend title.
  labs(
    title = "Genomic distribution of H3K9me3 and H4K20me3 peaks",
    subtitle = "Each bar represents the relative distribution of peak locations within one cell line and antibody",
    y = "Peaks (%)",
    fill = "Genomic location"
  )


# -----------------------
# Save outputs
# -----------------------

# Saves the absolute peak-count plot as a PDF.
ggsave(
  out_counts_pdf,
  plot = p_counts,
  width = 12,
  height = 9,
  useDingbats = FALSE
)


# Saves the same absolute peak-count plot as a high-resolution PNG.
ggsave(
  out_counts_png,
  plot = p_counts,
  width = 12,
  height = 9,
  dpi = 300
)


# Saves the percentage distribution plot as a PDF.
ggsave(
  out_percent_pdf,
  plot = p_percent,
  width = 12,
  height = 9,
  useDingbats = FALSE
)


# Saves the same percentage distribution plot as a high-resolution PNG.
ggsave(
  out_percent_png,
  plot = p_percent,
  width = 12,
  height = 9,
  dpi = 300
)


# Prints the paths of all generated tables and figures after the script finishes.
cat("\nFinished peak location plots.\n")
cat("Counts CSV:", out_counts_csv, "\n")
cat("Percent CSV:", out_percent_csv, "\n")
cat("Counts PDF:", out_counts_pdf, "\n")
cat("Counts PNG:", out_counts_png, "\n")
cat("Percent PDF:", out_percent_pdf, "\n")
cat("Percent PNG:", out_percent_png, "\n")