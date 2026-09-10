# This script annotates MACS2 ChIP-seq peaks according to their genomic location.
# ChIPseeker compares the peak coordinates with the mm10 gene annotation and assigns
# peaks to genomic features such as promoters, exons, introns, UTRs, downstream
# regions or intergenic regions.
#
# Each narrowPeak file, which contains the genomic regions identified as enriched by
# MACS2, is annotated separately. The complete annotation for every peak is saved for
# each sample, and the number of peaks assigned to each annotation category is collected
# in one summary table.
#
# Biologically, this is used to compare where H3K9me3 and H4K20me3 enrichment occurs
# in the genome across the different cell lines and genotypes. The annotation tables
# are also used later to create the genomic-location plots.


# Checks whether a CRAN package is available and installs it if needed.
# CRAN is the main repository used to distribute general R packages.
install_cran_if_missing <- function(pkg) {

  if (!requireNamespace(pkg, quietly = TRUE)) {

    message("Installing CRAN package: ", pkg)

    install.packages(
      pkg,
      repos = "https://cloud.r-project.org"
    )
  }
}


# Checks whether a Bioconductor package is available and installs it if needed.
install_bioc_if_missing <- function(pkg) {

  if (!requireNamespace(pkg, quietly = TRUE)) {

    message("Installing Bioconductor package: ", pkg)

    # Bioconductor provides R packages and reference data for bioinformatics
    # and biological data analysis. BiocManager is used to install and manage
    # these packages in R.
    if (!requireNamespace("BiocManager", quietly = TRUE)) {

      install.packages(
        "BiocManager",
        repos = "https://cloud.r-project.org"
      )
    }

    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )
  }
}


# Uses Conda to install R and system libraries required by some ChIPseeker dependencies.
# Conda manages software packages and their dependencies inside the active environment.
install_conda_packages <- function(pkgs) {

  # Finds the Conda program belonging to the active environment.
  conda <- Sys.which("conda")


  # Stops the script if Conda cannot be accessed from R.
  if (conda == "") {

    stop(
      "conda was not found from inside R. ",
      "Please run this script from an activated conda environment."
    )
  }


  # Prints which packages are being passed to Conda.
  message(
    "Installing/checking conda packages: ",
    paste(pkgs, collapse = " ")
  )


  # Runs the Conda installation command and stores its terminal output.
  status <- system2(
    conda,
    args = c(
      "install",
      "-y",
      "-c",
      "conda-forge",
      pkgs
    ),
    stdout = TRUE,
    stderr = TRUE
  )


  # Prints the Conda output so possible installation problems are visible.
  message(
    paste(status, collapse = "\n")
  )
}


# Checks that all packages needed for the peak annotation are available
# before the actual analysis starts.
ensure_required_packages <- function() {

  # ggplot2 provides plotting functions used by ChIPseeker and some of its dependencies.
  install_cran_if_missing("ggplot2")


  # These packages provide graphical, font and system libraries needed by
  # some of the R packages used together with ChIPseeker.
  conda_pkgs <- c(
    "r-ggiraph",
    "r-systemfonts",
    "r-gdtools",
    "cairo",
    "fontconfig",
    "freetype",
    "pkg-config"
  )


  # Uses Conda if one of the main graphical dependencies is still missing.
  if (
    !requireNamespace("ggiraph", quietly = TRUE) ||
    !requireNamespace("systemfonts", quietly = TRUE) ||
    !requireNamespace("gdtools", quietly = TRUE)
  ) {

    install_conda_packages(conda_pkgs)
  }


  # ChIPseeker assigns ChIP-seq peaks to genomic features and nearby genes.
  install_bioc_if_missing("ChIPseeker")


  # Provides a transcript database (TxDb) with the genomic positions of known
  # mouse genes and transcripts in the mm10 reference genome.
  install_bioc_if_missing(
    "TxDb.Mmusculus.UCSC.mm10.knownGene"
  )


  # Provides additional mouse gene information and gene identifiers,
  # which ChIPseeker can add to the peak annotation.
  install_bioc_if_missing("org.Mm.eg.db")


  # Lists the packages that must be available for the annotation itself.
  required <- c(
    "ChIPseeker",
    "TxDb.Mmusculus.UCSC.mm10.knownGene",
    "org.Mm.eg.db"
  )


  # Checks which required packages are still missing after the installation attempts.
  missing <- required[
    !vapply(
      required,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]


  # Stops the analysis if any required package is still unavailable.
  if (length(missing) > 0) {

    stop(
      "Still missing required packages after installation attempt: ",
      paste(missing, collapse = ", ")
    )
  }
}


# Runs the package check before the peak annotation starts.
ensure_required_packages()


# -----------------------
# Load libraries
# -----------------------

# Loads the packages used for peak annotation and mouse gene information.
# Startup messages are hidden so the console output stays easier to read.
suppressPackageStartupMessages({

  # Annotates ChIP-seq peaks relative to genomic features and nearby genes.
  library(ChIPseeker)

  # Provides the mm10 transcript database used as the genomic reference.
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)

  # Provides additional mouse gene annotation and identifiers.
  library(org.Mm.eg.db)
})


# =========================================================
# Settings
# =========================================================

# Uses the mm10 transcript database as the reference for assigning peaks
# to genomic features such as promoters, exons and introns.
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene


# Folder containing the MACS2 narrowPeak files.
peak_dir <- "03_macs2/peaks"


# Folder where the annotation tables are saved.
output_dir <- "05_analysis/04_annotation/annotation_tables"

dir.create(
  output_dir,
  showWarnings = FALSE
)


# Finds all MACS2 narrowPeak files in the input folder.
peak_files <- list.files(
  peak_dir,
  pattern = "_peaks\\.narrowPeak$",
  full.names = TRUE
)


# Stops the script if no peak files were found.
if (length(peak_files) == 0) {

  stop(
    "No *_peaks.narrowPeak files found in: ",
    peak_dir
  )
}


# Stores the annotation counts from each sample before they are combined.
summary_list <- list()


# =========================================================
# Annotate each peak file
# =========================================================

# Processes each MACS2 peak file separately.
for (peak_file in peak_files) {

  # Extracts the sample name from the peak filename.
  sample_name <- basename(peak_file)

  sample_name <- sub(
    "_peaks\\.narrowPeak$",
    "",
    sample_name
  )


  # Prints the sample that is currently being annotated.
  message(
    "Annotating: ",
    sample_name
  )


  # Assigns each peak to genomic features using the mm10 gene annotation.
  # The transcription start site (TSS) is the genomic position where transcription
  # of a gene begins. Peaks within 3 kb upstream or downstream of the TSS are
  # treated as promoter-associated for this annotation.
  peak_anno <- annotatePeak(
    peak_file,
    tssRegion = c(-3000, 3000),
    TxDb = txdb,
    annoDb = "org.Mm.eg.db"
  )


  # Converts the ChIPseeker result into a regular table with one row
  # containing the annotation information for each peak.
  annotation_df <- as.data.frame(
    peak_anno
  )


  # Saves the complete peak annotation table for this sample.
  write.csv(
    annotation_df,
    file = file.path(
      output_dir,
      paste0(
        sample_name,
        "_peak_annotation.csv"
      )
    ),
    row.names = FALSE
  )


  # Counts how many peaks were assigned to each annotation category.
  annotation_counts <- as.data.frame(
    table(annotation_df$annotation)
  )


  # Gives the two columns clearer names.
  colnames(annotation_counts) <- c(
    "annotation",
    "peak_count"
  )


  # Adds the sample name so counts from different samples can be identified later.
  annotation_counts$sample <- sample_name


  # Stores the annotation counts from this sample in the summary list.
  summary_list[[sample_name]] <- annotation_counts
}


# Combines the annotation counts from all samples into one data frame.
summary_df <- do.call(
  rbind,
  summary_list
)


# Keeps the sample name, annotation category and corresponding peak count.
summary_df <- summary_df[
  ,
  c(
    "sample",
    "annotation",
    "peak_count"
  )
]


# Saves the combined annotation summary for all samples.
write.csv(
  summary_df,
  file = file.path(
    output_dir,
    "all_peak_annotation_summary.csv"
  ),
  row.names = FALSE
)


message("Peak annotation finished.")
message("Results saved in: ", output_dir)