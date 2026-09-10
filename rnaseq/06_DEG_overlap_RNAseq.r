# This script compares differentially expressed genes (DEGs) across all knockout
# genotypes to identify expression changes that are shared between several knockouts
# or specific to individual genotypes.
#
# For each KO-versus-WT DESeq2 comparison, genes are classified as differentially
# expressed when they pass both the adjusted p-value and log2 fold-change thresholds.
# Separate DEG lists are created for all significant genes, upregulated genes and
# downregulated genes.
#
# The overlap between these gene sets is then summarized in two ways. A pairwise
# overlap matrix counts how many DEGs are shared between every pair of genotype
# comparisons, while UpSet plots show larger intersection patterns across multiple
# genotypes at the same time.
#
# This makes it possible to distinguish expression changes that occur repeatedly
# across several HP1 knockout backgrounds from changes that are more genotype-specific,
# and to compare these patterns with the H3K9 methyltransferase quadruple knockout.


# -----------------------
# Install missing packages
# -----------------------

# Defines the CRAN repository used to install general R packages.
# CRAN is the main online repository for R packages.
cran_repo <- "https://cloud.r-project.org"


# Checks whether ggplot2 is available and installs it from CRAN if needed.
if (!requireNamespace("ggplot2", quietly = TRUE)) {

  install.packages(
    "ggplot2",
    repos = cran_repo
  )
}


# Checks whether UpSetR is available and installs it from CRAN if needed.
if (!requireNamespace("UpSetR", quietly = TRUE)) {

  install.packages(
    "UpSetR",
    repos = cran_repo
  )
}


# -----------------------
# Load libraries
# -----------------------

# ggplot2 provides plotting functions and is also used internally by some
# visualization workflows in the analysis environment.
library(ggplot2)

# UpSetR is used to visualize intersections between DEG lists from
# multiple genotype comparisons.
library(UpSetR)


# -----------------------
# Settings
# -----------------------

# Defines the folder containing the DESeq2 comparison results and the folder
# where the DEG overlap analysis is saved.
in_dir <- "/home1/RNAseq/04_DESeq2_results"

out_dir <- "/home1/RNAseq/05_analysis/03_DEG_overlap"


# Creates the output folder if it does not already exist.
dir.create(
  out_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# Defines the thresholds used to classify genes as differentially expressed.
# Genes must have an adjusted p-value below 0.05 and an absolute log2 fold
# change greater than 1.
padj_cutoff <- 0.05

log2fc_cutoff <- 1


# Defines all knockout genotypes that are compared with WT.
comparisons <- c(
  "HP1a_KO",
  "HP1b_KO",
  "HP1g_KO",
  "HP1ab_DKO",
  "HP1ag_DKO",
  "HP1bg_DKO",
  "HP1abg_TKO",
  "H3K9MT_QKO"
)


# -----------------------
# Read DEG lists
# -----------------------

# Creates empty lists that will store all significant DEGs, upregulated DEGs
# and downregulated DEGs separately for each KO-versus-WT comparison.
deg_lists_all <- list()

deg_lists_up <- list()

deg_lists_down <- list()


# Reads the DESeq2 result file for each knockout and extracts the corresponding
# DEG lists according to the defined significance thresholds.
for (comparison in comparisons) {

  result_file <- file.path(
    in_dir,
    paste0(
      "DESeq2_",
      comparison,
      "_vs_WT_results.csv"
    )
  )


  # Skips the current comparison if its DESeq2 result file is missing and
  # prints a warning instead of stopping the complete analysis.
  if (!file.exists(result_file)) {

    warning(
      paste(
        "Result file not found and will be skipped:",
        result_file
      )
    )

    next
  }


  # Reads the DESeq2 result table for the current KO-versus-WT comparison.
  res <- read.csv(
    result_file,
    row.names = 1
  )


  # Removes genes without an adjusted p-value or log2 fold-change value.
  res <- res[
    !is.na(res$padj) &
    !is.na(res$log2FoldChange),
  ]


  # Extracts all significantly differentially expressed genes, independent
  # of whether expression is increased or decreased in the knockout.
  all_deg <- rownames(
    res[
      res$padj < padj_cutoff &
      abs(res$log2FoldChange) > log2fc_cutoff,
    ]
  )


  # Extracts genes that are significantly more highly expressed in the
  # knockout than in WT.
  up_deg <- rownames(
    res[
      res$padj < padj_cutoff &
      res$log2FoldChange > log2fc_cutoff,
    ]
  )


  # Extracts genes that are significantly less highly expressed in the
  # knockout than in WT.
  down_deg <- rownames(
    res[
      res$padj < padj_cutoff &
      res$log2FoldChange < -log2fc_cutoff,
    ]
  )


  # Stores the three DEG sets under the name of the current genotype comparison.
  deg_lists_all[[comparison]] <- all_deg

  deg_lists_up[[comparison]] <- up_deg

  deg_lists_down[[comparison]] <- down_deg
}


# Stops the script if none of the expected DESeq2 result files could be loaded.
if (length(deg_lists_all) == 0) {

  stop(
    paste(
      "No DESeq2 result files were found. Check file names in:",
      in_dir
    )
  )
}


# -----------------------
# Save DEG list sizes
# -----------------------

# Creates a summary table containing the number of all, upregulated and
# downregulated DEGs for every genotype comparison.
deg_size_table <- data.frame(
  comparison = names(deg_lists_all),
  all_DEGs = sapply(
    deg_lists_all,
    length
  ),
  upregulated = sapply(
    deg_lists_up,
    length
  ),
  downregulated = sapply(
    deg_lists_down,
    length
  )
)


# Saves the DEG counts used as input for the overlap analysis.
write.csv(
  deg_size_table,
  file = file.path(
    out_dir,
    "DEG_overlap_input_sizes.csv"
  ),
  row.names = FALSE
)


# -----------------------
# Save pairwise overlap matrix
# -----------------------

# Stores the names of all genotype comparisons that were successfully loaded.
all_names <- names(
  deg_lists_all
)


# Creates an empty square matrix in which rows and columns represent
# genotype comparisons and each cell will contain their number of shared DEGs.
overlap_matrix <- matrix(
  0,
  nrow = length(all_names),
  ncol = length(all_names)
)


# Labels the rows and columns with the genotype comparison names.
rownames(overlap_matrix) <- all_names

colnames(overlap_matrix) <- all_names


# Calculates the number of DEGs shared between every pair of genotype comparisons.
# The diagonal therefore contains the total DEG number of each comparison.
for (i in all_names) {

  for (j in all_names) {

    overlap_matrix[i, j] <- length(
      intersect(
        deg_lists_all[[i]],
        deg_lists_all[[j]]
      )
    )
  }
}


# Saves the pairwise DEG overlap counts as a CSV file.
write.csv(
  overlap_matrix,
  file = file.path(
    out_dir,
    "DEG_pairwise_overlap_matrix.csv"
  )
)


# -----------------------
# UpSet plot: all DEGs
# -----------------------

# Opens a PDF graphics device for the UpSet plot containing all significant DEGs.
pdf(
  file.path(
    out_dir,
    "DEG_overlap_UpSet_all_DEGs.pdf"
  ),
  width = 10,
  height = 6
)


# Creates an UpSet plot showing how many DEGs are unique to individual
# comparisons or shared across several knockout genotypes.
upset(
  fromList(
    deg_lists_all
  ),
  order.by = "freq",
  nsets = length(
    deg_lists_all
  ),
  mainbar.y.label = "DEG intersections",
  sets.x.label = "DEGs per comparison",
  main.bar.color = "#7EA6D9",
  sets.bar.color = "#B79AD6",
  matrix.color = "#55B7B1"
)


# Closes and saves the PDF graphics device.
dev.off()


# Opens a PNG graphics device for the same all-DEG UpSet plot.
png(
  file.path(
    out_dir,
    "DEG_overlap_UpSet_all_DEGs.png"
  ),
  width = 3000,
  height = 1800,
  res = 300
)


# Creates the same UpSet plot in high-resolution PNG format.
upset(
  fromList(
    deg_lists_all
  ),
  order.by = "freq",
  nsets = length(
    deg_lists_all
  ),
  mainbar.y.label = "DEG intersections",
  sets.x.label = "DEGs per comparison",
  main.bar.color = "#7EA6D9",
  sets.bar.color = "#B79AD6",
  matrix.color = "#55B7B1"
)


# Closes and saves the PNG graphics device.
dev.off()


# -----------------------
# UpSet plot: upregulated DEGs
# -----------------------

# Opens a PDF graphics device for the UpSet plot containing only genes
# that are significantly upregulated in the knockout compared with WT.
pdf(
  file.path(
    out_dir,
    "DEG_overlap_UpSet_upregulated_DEGs.pdf"
  ),
  width = 10,
  height = 6
)


# Shows which upregulated DEGs are unique to individual knockouts or shared
# between several genotype comparisons.
upset(
  fromList(
    deg_lists_up
  ),
  order.by = "freq",
  nsets = length(
    deg_lists_up
  ),
  mainbar.y.label = "Upregulated DEG intersections",
  sets.x.label = "Upregulated DEGs per comparison",
  main.bar.color = "#D65F8C",
  sets.bar.color = "#E8899C",
  matrix.color = "#D65F8C"
)


# Closes and saves the PDF graphics device.
dev.off()


# Opens a PNG graphics device for the same upregulated-DEG UpSet plot.
png(
  file.path(
    out_dir,
    "DEG_overlap_UpSet_upregulated_DEGs.png"
  ),
  width = 3000,
  height = 1800,
  res = 300
)


# Creates the same upregulated-DEG UpSet plot in high-resolution PNG format.
upset(
  fromList(
    deg_lists_up
  ),
  order.by = "freq",
  nsets = length(
    deg_lists_up
  ),
  mainbar.y.label = "Upregulated DEG intersections",
  sets.x.label = "Upregulated DEGs per comparison",
  main.bar.color = "#D65F8C",
  sets.bar.color = "#E8899C",
  matrix.color = "#D65F8C"
)


# Closes and saves the PNG graphics device.
dev.off()


# -----------------------
# UpSet plot: downregulated DEGs
# -----------------------

# Opens a PDF graphics device for the UpSet plot containing only genes
# that are significantly downregulated in the knockout compared with WT.
pdf(
  file.path(
    out_dir,
    "DEG_overlap_UpSet_downregulated_DEGs.pdf"
  ),
  width = 10,
  height = 6
)


# Shows which downregulated DEGs are unique to individual knockouts or shared
# between several genotype comparisons.
upset(
  fromList(
    deg_lists_down
  ),
  order.by = "freq",
  nsets = length(
    deg_lists_down
  ),
  mainbar.y.label = "Downregulated DEG intersections",
  sets.x.label = "Downregulated DEGs per comparison",
  main.bar.color = "#55B7B1",
  sets.bar.color = "#BFE9E4",
  matrix.color = "#55B7B1"
)


# Closes and saves the PDF graphics device.
dev.off()


# Opens a PNG graphics device for the same downregulated-DEG UpSet plot.
png(
  file.path(
    out_dir,
    "DEG_overlap_UpSet_downregulated_DEGs.png"
  ),
  width = 3000,
  height = 1800,
  res = 300
)


# Creates the same downregulated-DEG UpSet plot in high-resolution PNG format.
upset(
  fromList(
    deg_lists_down
  ),
  order.by = "freq",
  nsets = length(
    deg_lists_down
  ),
  mainbar.y.label = "Downregulated DEG intersections",
  sets.x.label = "Downregulated DEGs per comparison",
  main.bar.color = "#55B7B1",
  sets.bar.color = "#BFE9E4",
  matrix.color = "#55B7B1"
)


# Closes and saves the PNG graphics device.
dev.off()


# Confirms that the DEG overlap analysis was completed and prints the output location.
cat(
  "\nDEG overlap analysis finished.\n"
)

cat(
  "Results saved in:",
  out_dir,
  "\n"
)