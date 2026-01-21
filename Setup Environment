#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# Environment Setup Script for Microwave Stimulation Analysis Suite
# -----------------------------------------------------------------------------
# This script sets up the complete R environment required for all analyses
# Run this script first before using any analysis modules
# -----------------------------------------------------------------------------

cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║      Microwave Stimulation Analysis Suite - Environment Setup        ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

# Record start time
start_time <- Sys.time()
cat("Setup started at:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n\n")

# =============================================================================
# SECTION 1: CHECK SYSTEM REQUIREMENTS
# =============================================================================

cat("1. Checking system requirements...\n")

# Check R version
r_version <- getRversion()
required_r_version <- "4.1.0"
cat("   R version:", as.character(r_version), "\n")

if (r_version < required_r_version) {
  cat("   ⚠️  WARNING: R version should be at least", required_r_version, "\n")
  cat("   Consider updating R from https://cloud.r-project.org/\n")
} else {
  cat("   ✓ R version OK\n")
}

# Check available memory
if (.Platform$OS.type == "windows") {
  mem_info <- memory.size(max = FALSE)
  cat("   Available memory:", round(mem_info / 1024^2, 1), "MB\n")
} else {
  # Unix-like systems
  mem_info <- system("grep MemTotal /proc/meminfo", intern = TRUE)
  if (length(mem_info) > 0) {
    mem_mb <- as.numeric(gsub("[^0-9]", "", mem_info)) / 1024
    cat("   Total memory:", round(mem_mb, 1), "MB\n")
  }
}

# Check working directory
cat("   Working directory:", getwd(), "\n\n")

# =============================================================================
# SECTION 2: CRAN PACKAGES
# =============================================================================

cat("2. Installing/updating CRAN packages...\n")

cran_packages <- c(
  # Core data manipulation
  "tidyverse", "dplyr", "tidyr", "readr", "purrr", "stringr", "forcats",
  "data.table", "vroom", "readxl", "writexl", "openxlsx",
  
  # Visualization
  "ggplot2", "ggpubr", "cowplot", "patchwork", "ggrepel", "ggforce",
  "ggdist", "ggExtra", "ggthemes", "RColorBrewer", "viridis", "wesanderson",
  "scales", "gridExtra", "grid", "corrplot", "pheatmap", "heatmaply",
  
  # Statistical analysis
  "rstatix", "broom", "emmeans", "multcomp", "car", "MASS", "lme4",
  "nlme", "survival", "survminer", "coin", "boot", "simpleboot",
  "moments", "fitdistrplus", "goftest",
  
  # Advanced analysis
  "factoextra", "cluster", "NbClust", "dbscan", "umap", "tsne", "Rtsne",
  "seriation", "ape", "dendextend", "phytools", "dynamicTreeCut",
  "randomForest", "caret", "glmnet", "xgboost", "ranger", "e1071",
  "kernlab", "naivebayes",
  
  # Time series/signal processing
  "zoo", "xts", "signal", "pracma", "changepoint", "strucchange",
  
  # Utilities
  "here", "fs", "glue", "magrittr", "knitr", "rmarkdown", "kableExtra",
  "bookdown", "tinytex", "DT", "shiny", "plotly", "htmlwidgets",
  "webshot", "webshot2", "reticulate", "devtools", "remotes", "pak",
  "sessioninfo", "renv", "packrat", "checkpoint",
  
  # Specialized packages
  "ggbeeswarm", "ggridges", "ggalluvial", "ggstream", "gganimate",
  "ggmap", "sf", "maps", "mapdata",
  "ComplexHeatmap", "circlize", "UpSetR", "VennDiagram", "eulerr",
  "pROC", "ROCR", "PRROC",
  "missForest", "mice", "VIM", "naniar",
  "pwr", "MBESS", "effectsize",
  "modelr", "rsample", "yardstick", "dials", "tune", "workflows",
  "recipes", "embed", "textrecipes", "themis", "butcher", "probably"
)

# Function to install packages with error handling
install_with_retry <- function(packages, repos = "https://cloud.r-project.org/") {
  installed <- character()
  failed <- character()
  
  for (pkg in packages) {
    cat("   Installing", pkg, "... ")
    
    # Check if already installed
    if (!requireNamespace(pkg, quietly = TRUE)) {
      # Try to install
      result <- tryCatch({
        install.packages(pkg, repos = repos, dependencies = TRUE, quiet = TRUE)
        TRUE
      }, error = function(e) FALSE)
      
      if (result) {
        cat("✓\n")
        installed <- c(installed, pkg)
      } else {
        cat("✗\n")
        failed <- c(failed, pkg)
      }
    } else {
      cat("already installed ✓\n")
    }
  }
  
  return(list(installed = installed, failed = failed))
}

# Install CRAN packages
cran_result <- install_with_retry(cran_packages)

cat("\n   CRAN packages summary:\n")
cat("   - Successfully installed:", length(cran_result$installed), "\n")
cat("   - Failed to install:", length(cran_result$failed), "\n")
if (length(cran_result$failed) > 0) {
  cat("   Failed packages:", paste(cran_result$failed, collapse = ", "), "\n")
}

# =============================================================================
# SECTION 3: BIOCONDUCTOR PACKAGES
# =============================================================================

cat("\n3. Installing/updating Bioconductor packages...\n")

# Install BiocManager if not present
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("   Installing BiocManager... ")
  install.packages("BiocManager", repos = "https://cloud.r-project.org/")
  cat("✓\n")
}

# Bioconductor packages
bioc_packages <- c(
  # RNA-seq analysis
  "limma", "edgeR", "DESeq2", "DEXSeq", "tximport", "tximeta",
  "sva", "RUVSeq", "clusterProfiler", "enrichplot", "DOSE",
  "gage", "gageData", "pathview", "reactome.db", "reactomePA",
  "org.Hs.eg.db", "org.Mm.eg.db",
  "fgsea", "msigdbr", "GSVA", "GSEABase", "GO.db", "KEGG.db",
  "BiocParallel", "BiocFileCache", "AnnotationDbi", "ensembldb",
  "GenomicFeatures", "GenomicAlignments", "Rsamtools", "Rsubread",
  
  # Single cell/specialized
  "Seurat", "SingleCellExperiment", "scater", "scran", "DropletUtils",
  "celda", "zinbwave", "MAST", "monocle", "slingshot", "tradeSeq",
  "BUSpaRse", "DropletQC",
  
  # Visualization
  "ComplexHeatmap", "ggbio", "Gviz", "karyoploteR", "trackViewer",
  "CHETAH", "scRepertoire", "scCustomize",
  
  # Utilities
  "BiocCheck", "BiocStyle", "BiocGenerics", "Biobase", "SummarizedExperiment",
  "MultiAssayExperiment", "ExperimentHub", "AnnotationHub", "ensembldb",
  "BSgenome", "rtracklayer", "GenomicRanges", "IRanges", "VariantAnnotation"
)

# Function to install Bioconductor packages
install_bioc <- function(packages) {
  installed <- character()
  failed <- character()
  
  for (pkg in packages) {
    cat("   Installing", pkg, "... ")
    
    if (!requireNamespace(pkg, quietly = TRUE)) {
      result <- tryCatch({
        BiocManager::install(pkg, update = FALSE, ask = FALSE, quiet = TRUE)
        TRUE
      }, error = function(e) FALSE)
      
      if (result) {
        cat("✓\n")
        installed <- c(installed, pkg)
      } else {
        cat("✗\n")
        failed <- c(failed, pkg)
      }
    } else {
      cat("already installed ✓\n")
    }
  }
  
  return(list(installed = installed, failed = failed))
}

# Install Bioconductor packages
bioc_result <- install_bioc(bioc_packages)

cat("\n   Bioconductor packages summary:\n")
cat("   - Successfully installed:", length(bioc_result$installed), "\n")
cat("   - Failed to install:", length(bioc_result$failed), "\n")

# =============================================================================
# SECTION 4: GITHUB PACKAGES
# =============================================================================

cat("\n4. Installing GitHub packages...\n")

github_packages <- list(
  # Visualization
  c("smin95/smplot2", "smplot2"),
  c("nsgrantham/ggdark", "ggdark"),
  c("hrbrmstr/streamgraph", "streamgraph"),
  c("davidsjoberg/ggstream", "ggstream"),
  c("exaexa/scattermore", "scattermore"),
  c("cardiomoon/moonBook", "moonBook"),
  c("cardiomoon/webr", "webr"),
  c("kassambara/ggcorrplot", "ggcorrplot"),
  c("kassambara/ggsurvfit", "ggsurvfit"),
  c("const-ae/ggsignif", "ggsignif"),
  c("rstudio/gt", "gt"),
  
  # Analysis utilities
  c("ropensci/plotly", "plotly"),
  c("rstudio/reticulate", "reticulate"),
  c("r-lib/testthat", "testthat"),
  c("r-lib/usethis", "usethis"),
  c("r-lib/pkgdown", "pkgdown"),
  c("r-lib/styler", "styler"),
  c("r-lib/lintr", "lintr"),
  c("ThinkR-open/attachment", "attachment"),
  
  # Specialized
  c("satijalab/seurat", "Seurat"),
  c("GreenleafLab/ArchR", "ArchR"),
  c("immunogenomics/harmony", "harmony")
)

# Function to install GitHub packages
install_github_pkgs <- function(package_list) {
  installed <- character()
  failed <- character()
  
  for (pkg_info in package_list) {
    repo <- pkg_info[1]
    pkg_name <- pkg_info[2]
    
    cat("   Installing", pkg_name, "from", repo, "... ")
    
    if (!requireNamespace(pkg_name, quietly = TRUE)) {
      result <- tryCatch({
        remotes::install_github(repo, upgrade = "never", quiet = TRUE)
        TRUE
      }, error = function(e) FALSE)
      
      if (result) {
        cat("✓\n")
        installed <- c(installed, pkg_name)
      } else {
        cat("✗\n")
        failed <- c(failed, pkg_name)
      }
    } else {
      cat("already installed ✓\n")
    }
  }
  
  return(list(installed = installed, failed = failed))
}

# Install GitHub packages
github_result <- install_github_pkgs(github_packages)

cat("\n   GitHub packages summary:\n")
cat("   - Successfully installed:", length(github_result$installed), "\n")
cat("   - Failed to install:", length(github_result$failed), "\n")

# =============================================================================
# SECTION 5: VERIFICATION AND TESTING
# =============================================================================

cat("\n5. Verifying installation...\n")

# Test loading core packages
core_test_packages <- c(
  "tidyverse", "ggplot2", "dplyr", "tidyr",
  "limma", "edgeR", "clusterProfiler",
  "ggpubr", "patchwork", "rstatix"
)

loading_issues <- character()
for (pkg in core_test_packages) {
  cat("   Testing", pkg, "... ")
  success <- tryCatch({
    library(pkg, character.only = TRUE, quietly = TRUE)
    TRUE
  }, error = function(e) FALSE)
  
  if (success) {
    cat("✓\n")
  } else {
    cat("✗\n")
    loading_issues <- c(loading_issues, pkg)
  }
}

# Create requirements file
cat("\n6. Creating requirements file...\n")

# Save package versions
package_versions <- data.frame(
  Package = character(),
  Version = character(),
  Source = character(),
  stringsAsFactors = FALSE
)

all_packages <- c(cran_packages, bioc_packages, sapply(github_packages, `[`, 2))

for (pkg in all_packages) {
  version <- tryCatch({
    as.character(packageVersion(pkg))
  }, error = function(e) "Not installed")
  
  source <- if (pkg %in% cran_packages) {
    "CRAN"
  } else if (pkg %in% bioc_packages) {
    "Bioconductor"
  } else {
    "GitHub"
  }
  
  package_versions <- rbind(package_versions, 
                            data.frame(Package = pkg, 
                                       Version = version, 
                                       Source = source))
}

write.csv(package_versions, "package_versions.csv", row.names = FALSE)
cat("   ✓ Package versions saved to: package_versions.csv\n")

# Create renv lockfile (optional)
if (requireNamespace("renv", quietly = TRUE)) {
  cat("   Creating renv lockfile... ")
  renv::snapshot(lockfile = "renv.lock", prompt = FALSE)
  cat("✓\n")
}

# =============================================================================
# SECTION 6: CREATE PROJECT STRUCTURE
# =============================================================================

cat("\n7. Setting up project structure...\n")

# Define directory structure
dirs_to_create <- c(
  "data",
  "data/morphology",
  "data/calcium",
  "data/cytotoxicity", 
  "data/transcriptomics",
  "data/example",
  "results",
  "results/morphology",
  "results/calcium",
  "results/cytotoxicity",
  "results/transcriptomics",
  "reports",
  "figures",
  "figures/publication",
  "figures/exploratory",
  "scripts",
  "scripts/utilities",
  "scripts/analysis",
  "config",
  "logs",
  "temp"
)

# Create directories
created_dirs <- character()
for (dir in dirs_to_create) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    created_dirs <- c(created_dirs, dir)
  }
}

cat("   Created", length(created_dirs), "directories\n")

# Create README files
readme_content <- "# Analysis Directory\n\nThis directory contains analysis outputs.\n\n## Structure\n- `data/`: Raw and processed data\n- `results/`: Analysis results\n- `figures/`: Generated plots\n- `reports/`: HTML/PDF reports\n- `scripts/`: R scripts\n"

writeLines(readme_content, "README_STRUCTURE.md")
cat("   ✓ Project structure documentation created\n")

# =============================================================================
# SECTION 7: CREATE CONFIGURATION FILES
# =============================================================================

cat("\n8. Creating configuration files...\n")

# Create .Rprofile
rprofile_content <- '# Microwave Stimulation Analysis Suite - .Rprofile
# This file runs when R starts in this project

cat("\\n=============================================\\n")
cat("  Microwave Stimulation Analysis Suite\\n")
cat("=============================================\\n\\n")

# Set project options
options(
  stringsAsFactors = FALSE,
  scipen = 999,           # Disable scientific notation
  digits = 4,             # Number of digits to display
  warn = -1,              # Turn off warnings (use with caution)
  error = NULL,
  repos = "https://cloud.r-project.org/"
)

# Load commonly used packages on startup
if (interactive()) {
  suppressPackageStartupMessages({
    library(here)
    library(tidyverse)
    library(ggplot2)
  })
  
  cat("Loaded: here, tidyverse, ggplot2\\n")
}

# Set ggplot2 theme
theme_set(theme_minimal(base_size = 11) +
            theme(plot.title = element_text(face = "bold", hjust = 0.5),
                  legend.position = "right",
                  panel.grid.minor = element_blank()))

# Custom functions
cat("Type `setup_help()` for available helper functions\\n")

setup_help <- function() {
  cat("\\nAvailable helper functions:\\n")
  cat("  setup_help()          - Show this help\\n")
  cat("  check_env()           - Check environment setup\\n")
  cat("  update_packages()     - Update all packages\\n")
  cat("  clear_temp()          - Clear temporary files\\n")
  cat("\\nAnalysis modules:\\n")
  cat("  run_morphology()      - Run morphological analysis\\n")
  cat("  run_calcium()         - Run calcium imaging analysis\\n")
  cat("  run_cytotoxicity()    - Run LDH assay analysis\\n")
  cat("  run_transcriptomics() - Run RNA-seq analysis\\n")
}

check_env <- function() {
  cat("\\nEnvironment check:\\n")
  cat("  R version:", as.character(getRversion()), "\\n")
  cat("  Platform:", R.version$platform, "\\n")
  cat("  Working directory:", getwd(), "\\n")
  cat("  Project root:", here::here(), "\\n")
}

update_packages <- function() {
  cat("Updating packages...\\n")
  update.packages(ask = FALSE, checkBuilt = TRUE)
  if (requireNamespace("BiocManager", quietly = TRUE)) {
    BiocManager::install(ask = FALSE)
  }
  cat("Done!\\n")
}

clear_temp <- function() {
  temp_dir <- "temp"
  if (dir.exists(temp_dir)) {
    unlink(temp_dir, recursive = TRUE)
    dir.create(temp_dir)
    cat("Temporary directory cleared\\n")
  }
}

# Project-specific paths
paths <- list(
  data = "data",
  results = "results",
  figures = "figures",
  reports = "reports",
  scripts = "scripts"
)

# Make paths available globally
list2env(paths, envir = .GlobalEnv)

cat("\\nProject ready! Type `setup_help()` for available commands.\\n\\n")
'

writeLines(rprofile_content, ".Rprofile")
cat("   ✓ Created .Rprofile\n")

# Create .gitignore
gitignore_content <- '# R
.Rhistory
.RData
.Ruserdata
.Rproj.user/
*.Rproj

# Temporary files
*.log
*.aux
*.out
*.toc
*.lot
*.lof
*.blg
*.bbl
*.synctex.gz
*.fdb_latexmk
*.fls

# Data files
*.csv
*.xlsx
*.xls
*.rds
*.rda
*.Rdata

# Results
results/*
!results/README.md
figures/*
!figures/README.md

# Environment
.Renviron
renv/
.Rprofile
.Rproj.user
*.Rproj

# OS files
.DS_Store
Thumbs.db
*.lnk

# IDE
.idea/
.vscode/
*.swp
*.swo

# Logs
logs/*
!logs/README.md

# Temporary
temp/*
!temp/README.md

# Large files
*.pdf
*.html
*.zip
*.tar.gz

# Package versions
package_versions.csv
'

writeLines(gitignore_content, ".gitignore")
cat("   ✓ Created .gitignore\n")

# Create configuration file
config_content <- '# Microwave Stimulation Analysis Suite - Configuration
# This file contains project-wide configuration settings

# ============================================================================
# GENERAL SETTINGS
# ============================================================================

# Project information
project_name <- "Microwave Stimulation Analysis"
project_version <- "1.0.0"
author <- "Your Name"
institution <- "Your Institution"

# Analysis settings
seed <- 42  # For reproducibility
cores <- parallel::detectCores() - 1  # Leave one core free

# ============================================================================
# VISUALIZATION SETTINGS
# ============================================================================

# Color palettes
palettes <- list(
  
  # Treatment colors
  treatment = c(
    "OFF" = "darkgrey",
    "ON" = "#E50000",
    "Sham" = "darkgrey",
    "Stimulation" = "#E50000",
    "Control" = "#1F77B4",
    "Heat" = "#FF7F0E"
  ),
  
  # Cell type colors
  cell_type = c(
    "Astrocyte" = "#016FB9",
    "Microglia" = "#B6C649",
    "GBM" = "#FB5607"
  ),
  
  # Region colors
  region = c(
    "Center" = "#FB5607",
    "Periphery" = "#016FB9"
  ),
  
  # Recovery time colors
  recovery = c(
    "0" = "#FB5607",
    "24" = "#B6C649",
    "48" = "#016FB9"
  ),
  
  # Volcano plot colors
  volcano = c(
    "Up-regulated" = "#FD841F",
    "Down-regulated" = "#38E54D",
    "Stable" = "grey90"
  )
)

# Plotting theme
plot_theme <- function(base_size = 11, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold")
    )
}

# Figure dimensions
figure_sizes <- list(
  small = c(width = 3.5, height = 3),
  medium = c(width = 5, height = 4),
  large = c(width = 7, height = 5),
  wide = c(width = 10, height = 6),
  publication = c(width = 8.5, height = 11)  # Letter size
)

# ============================================================================
# ANALYSIS PARAMETERS
# ============================================================================

# Statistical thresholds
thresholds <- list(
  p_value = 0.05,
  fdr = 0.05,
  fold_change = 1.0,
  min_cpm = 1.0,
  min_samples = 2
)

# Calcium imaging parameters
calcium <- list(
  recording_time = 3,           # minutes
  cutoff_fixed = 0.2,           # activity threshold
  delta_multiplier = 2,         # peak detection sensitivity
  window_size = 50,             # baseline correction window
  skip_correction = TRUE        # skip additional correction
)

# Morphology analysis parameters
morphology <- list(
  company_colors = c("#E50000", "#008A8A", "#AF0076", "#E56800", 
                    "#1717A0", "#E5AC00", "#00B700"),
  optimal_clusters = 6,
  pval_threshold = 0.05
)

# ============================================================================
# FILE PATHS
# ============================================================================

# Base directories (relative to project root)
paths <- list(
  data = "data",
  results = "results",
  figures = "figures",
  reports = "reports",
  scripts = "scripts",
  logs = "logs",
  temp = "temp"
)

# Data subdirectories
paths$data_subdirs <- list(
  morphology = file.path(paths$data, "morphology"),
  calcium = file.path(paths$data, "calcium"),
  cytotoxicity = file.path(paths$data, "cytotoxicity"),
  transcriptomics = file.path(paths$data, "transcriptomics"),
  example = file.path(paths$data, "example")
)

# Results subdirectories
paths$results_subdirs <- list(
  morphology = file.path(paths$results, "morphology"),
  calcium = file.path(paths$results, "calcium"),
  cytotoxicity = file.path(paths$results, "cytotoxicity"),
  transcriptomics = file.path(paths$results, "transcriptomics")
)

# ============================================================================
# EXPORT SETTINGS
# ============================================================================

export <- list(
  formats = c("pdf", "png", "svg"),
  dpi = 300,
  units = "in",
  device = "pdf"
)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# Function to get full path
get_path <- function(type, subdir = NULL, filename = NULL) {
  if (!is.null(subdir) && subdir %in% names(paths$data_subdirs)) {
    base_path <- paths$data_subdirs[[subdir]]
  } else if (!is.null(subdir) && subdir %in% names(paths$results_subdirs)) {
    base_path <- paths$results_subdirs[[subdir]]
  } else {
    base_path <- paths[[type]]
  }
  
  if (!is.null(filename)) {
    file.path(base_path, filename)
  } else {
    base_path
  }
}

# Function to save plot with consistent settings
save_plot <- function(plot, filename, size = "medium", 
                      subdir = NULL, format = "pdf") {
  
  if (!size %in% names(figure_sizes)) {
    size <- "medium"
  }
  
  dims <- figure_sizes[[size]]
  
  if (!is.null(subdir) && dir.exists(subdir)) {
    filename <- file.path(subdir, filename)
  }
  
  full_filename <- paste0(filename, ".", format)
  
  ggsave(full_filename, plot, 
         width = dims["width"], 
         height = dims["height"],
         dpi = export$dpi,
         units = export$units,
         device = format)
  
  cat("Plot saved to:", full_filename, "\\n")
}

# ============================================================================
# PRINT CONFIGURATION SUMMARY
# ============================================================================

cat_config_summary <- function() {
  cat("\\nConfiguration Summary:\\n")
  cat("=====================\\n")
  cat("Project:", project_name, "v", project_version, "\\n")
  cat("Author:", author, "\\n")
  cat("Institution:", institution, "\\n")
  cat("Seed:", seed, "\\n")
  cat("Available cores:", cores, "\\n")
  cat("\\nPaths:\\n")
  for (name in names(paths)) {
    if (is.character(paths[[name]])) {
      cat("  ", name, ":", paths[[name]], "\\n")
    }
  }
}

# Print summary when sourced
if (sys.nframe() == 0) {
  cat_config_summary()
}
'

writeLines(config_content, "config/settings.R")
cat("   ✓ Created config/settings.R\n")

# =============================================================================
# SECTION 8: CREATE EXAMPLE DATA AND SCRIPTS
# =============================================================================

cat("\n9. Creating example scripts...\n")

# Create example analysis script
example_script <- '# Microwave Stimulation Analysis - Example Script
# This script demonstrates how to use the analysis pipeline

# Load configuration
source("config/settings.R")

# Set random seed for reproducibility
set.seed(seed)

# Load required libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)

# -----------------------------------------------------------------------------
# EXAMPLE 1: MORPHOLOGICAL ANALYSIS
# -----------------------------------------------------------------------------

run_morphology_example <- function() {
  cat("Running morphological analysis example...\\n")
  
  # Create example data
  set.seed(42)
  example_data <- tibble(
    Condition = rep(c("G0", "G1", "G2", "G3"), each = 50),
    Area = c(rnorm(50, 1000, 200), rnorm(50, 1200, 250),
             rnorm(50, 900, 180), rnorm(50, 1100, 220)),
    Perimeter = c(rnorm(50, 150, 30), rnorm(50, 180, 35),
                  rnorm(50, 140, 28), rnorm(50, 170, 32)),
    Circularity = c(rnorm(50, 0.7, 0.1), rnorm(50, 0.6, 0.12),
                    rnorm(50, 0.75, 0.09), rnorm(50, 0.65, 0.11))
  )
  
  # Add metadata
  example_data <- example_data %>%
    mutate(
      MW = case_when(
        Condition %in% c("G0", "G2") ~ "ON",
        Condition %in% c("G1", "G3") ~ "OFF",
        TRUE ~ "Unknown"
      ),
      Expo_Time = case_when(
        Condition %in% c("G0", "G1") ~ "5",
        Condition %in% c("G2", "G3") ~ "30",
        TRUE ~ "Unknown"
      )
    )
  
  # Create plot
  p <- ggplot(example_data, aes(x = MW, y = Area, fill = MW)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.5, size = 1) +
    facet_wrap(~ Expo_Time) +
    scale_fill_manual(values = palettes$treatment) +
    labs(
      title = "Example: Cell Area by Condition",
      x = "Microwave Status",
      y = "Cell Area (pixels)"
    ) +
    plot_theme() +
    stat_compare_means(method = "t.test", label = "p.format")
  
  # Save plot
  save_plot(p, "example_morphology", size = "medium")
  
  cat("Morphology example complete! Check figures/ directory.\\n")
}

# -----------------------------------------------------------------------------
# EXAMPLE 2: CALCIUM IMAGING ANALYSIS
# -----------------------------------------------------------------------------

run_calcium_example <- function() {
  cat("Running calcium imaging example...\\n")
  
  # Create example calcium trace
  set.seed(42)
  time_points <- 1:300
  baseline <- 100
  events <- c(50, 150, 250)
  
  # Simulate calcium trace with events
  trace <- baseline + rnorm(length(time_points), 0, 5)
  
  for (event in events) {
    event_window <- event:(event + 20)
    trace[event_window] <- trace[event_window] + 
      dnorm(seq(-2, 2, length.out = 21)) * 50
  }
  
  # Create data frame
  calcium_data <- tibble(
    Time = time_points,
    Fluorescence = trace,
    Condition = "Example"
  )
  
  # Create plot
  p <- ggplot(calcium_data, aes(x = Time, y = Fluorescence)) +
    geom_line(color = "#016FB9", size = 1) +
    geom_vline(xintercept = events, linetype = "dashed", 
               color = "#E50000", alpha = 0.5) +
    labs(
      title = "Example: Calcium Trace",
      x = "Time (frames)",
      y = "Fluorescence (AU)"
    ) +
    plot_theme()
  
  # Save plot
  save_plot(p, "example_calcium", size = "wide")
  
  cat("Calcium imaging example complete!\\n")
}

# -----------------------------------------------------------------------------
# MAIN EXECUTION
# -----------------------------------------------------------------------------

main <- function() {
  cat("\\n=============================================\\n")
  cat("  EXAMPLE ANALYSIS SCRIPT\\n")
  cat("=============================================\\n\\n")
  
  # Create necessary directories
  dir.create("figures/example", showWarnings = FALSE, recursive = TRUE)
  
  # Run examples
  run_morphology_example()
  run_calcium_example()
  
  cat("\\nAll examples completed successfully!\\n")
  cat("Check the figures/example directory for outputs.\\n")
}

# Run if script is executed directly
if (sys.nframe() == 0) {
  main()
}
'

writeLines(example_script, "scripts/example_analysis.R")
cat("   ✓ Created scripts/example_analysis.R\n")

# Create utility script
utility_script <- '# Microwave Stimulation Analysis - Utility Functions
# This script contains helper functions for the analysis pipeline

# -----------------------------------------------------------------------------
# DATA VALIDATION FUNCTIONS
# -----------------------------------------------------------------------------

# Check if required columns exist
check_columns <- function(data, required_cols) {
  missing <- setdiff(required_cols, colnames(data))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

# Validate experimental design
validate_design <- function(metadata, condition_col = "Condition") {
  if (!condition_col %in% colnames(metadata)) {
    stop("Condition column not found: ", condition_col)
  }
  
  conditions <- unique(metadata[[condition_col]])
  cat("Found", length(conditions), "conditions:", 
      paste(conditions, collapse = ", "), "\\n")
  
  # Check for replicates
  rep_counts <- metadata %>%
    group_by(!!sym(condition_col)) %>%
    summarise(n = n(), .groups = "drop")
  
  if (any(rep_counts$n < 2)) {
    warning("Some conditions have less than 2 replicates")
  }
  
  return(rep_counts)
}

# -----------------------------------------------------------------------------
# DATA TRANSFORMATION FUNCTIONS
# -----------------------------------------------------------------------------

# Normalize to control
normalize_to_control <- function(data, value_col, group_col, 
                                 control_value = NULL) {
  if (is.null(control_value)) {
    control_mean <- data %>%
      filter(!!sym(group_col) == "Control") %>%
      pull(!!sym(value_col)) %>%
      mean(na.rm = TRUE)
  } else {
    control_mean <- control_value
  }
  
  data <- data %>%
    mutate(!!paste0("norm_", value_col) := !!sym(value_col) / control_mean)
  
  return(data)
}

# Calculate z-scores
calculate_zscore <- function(data, value_col, group_col = NULL) {
  if (is.null(group_col)) {
    data <- data %>%
      mutate(!!paste0("z_", value_col) := 
               (!!sym(value_col) - mean(!!sym(value_col), na.rm = TRUE)) / 
               sd(!!sym(value_col), na.rm = TRUE))
  } else {
    data <- data %>%
      group_by(!!sym(group_col)) %>%
      mutate(!!paste0("z_", value_col) := 
               (!!sym(value_col) - mean(!!sym(value_col), na.rm = TRUE)) / 
               sd(!!sym(value_col), na.rm = TRUE)) %>%
      ungroup()
  }
  
  return(data)
}

# -----------------------------------------------------------------------------
# STATISTICAL HELPER FUNCTIONS
# -----------------------------------------------------------------------------

# Perform multiple testing correction
correct_multiple_testing <- function(p_values, method = "BH") {
  p_adj <- p.adjust(p_values, method = method)
  return(p_adj)
}

# Calculate effect size (Cohen\'s d)
calculate_cohens_d <- function(x, y) {
  n1 <- length(x)
  n2 <- length(y)
  m1 <- mean(x, na.rm = TRUE)
  m2 <- mean(y, na.rm = TRUE)
  s1 <- var(x, na.rm = TRUE)
  s2 <- var(y, na.rm = TRUE)
  
  s_pooled <- sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
  d <- (m1 - m2) / s_pooled
  
  return(d)
}

# -----------------------------------------------------------------------------
# PLOTTING HELPER FUNCTIONS
# -----------------------------------------------------------------------------

# Create publication-ready theme
theme_publication <- function(base_size = 11, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = base_size + 2),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40", size = base_size),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(size = base_size - 1),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size - 1),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold", size = base_size),
      panel.grid.major = element_line(color = "grey90", size = 0.2),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey70", size = 0.5),
      plot.margin = unit(c(1, 1, 1, 1), "cm")
    )
}

# Add significance stars to plot
add_significance_stars <- function(p_value) {
  stars <- case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ "ns"
  )
  return(stars)
}

# -----------------------------------------------------------------------------
# FILE MANAGEMENT FUNCTIONS
# -----------------------------------------------------------------------------

# Create timestamped filename
timestamp_filename <- function(base_name, extension, 
                               timestamp = NULL, 
                               include_time = TRUE) {
  if (is.null(timestamp)) {
    timestamp <- Sys.time()
  }
  
  if (include_time) {
    time_str <- format(timestamp, "%Y%m%d_%H%M%S")
  } else {
    time_str <- format(timestamp, "%Y%m%d")
  }
  
  filename <- paste0(base_name, "_", time_str, ".", extension)
  return(filename)
}

# Save analysis results
save_results <- function(results, name, directory = "results") {
  # Create directory if it doesn\'t exist
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
  }
  
  # Generate filename
  filename <- timestamp_filename(name, "rds")
  full_path <- file.path(directory, filename)
  
  # Save results
  saveRDS(results, full_path)
  
  # Also save as CSV if it\'s a data frame
  if (is.data.frame(results)) {
    csv_filename <- timestamp_filename(name, "csv")
    csv_path <- file.path(directory, csv_filename)
    write.csv(results, csv_path, row.names = FALSE)
  }
  
  cat("Results saved to:", full_path, "\\n")
  return(full_path)
}

# -----------------------------------------------------------------------------
# LOGGING FUNCTIONS
# -----------------------------------------------------------------------------

# Create log file
setup_logging <- function(log_dir = "logs", session_id = NULL) {
  if (!dir.exists(log_dir)) {
    dir.create(log_dir, recursive = TRUE)
  }
  
  if (is.null(session_id)) {
    session_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
  }
  
  log_file <- file.path(log_dir, paste0("analysis_", session_id, ".log"))
  
  # Create log connection
  log_con <- file(log_file, open = "a")
  sink(log_con, split = TRUE)
  
  cat("\\n=============================================\\n")
  cat("Analysis Log - Session:", session_id, "\\n")
  cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n")
  cat("=============================================\\n\\n")
  
  return(list(con = log_con, file = log_file, id = session_id))
}

# Close logging
close_logging <- function(log_session) {
  cat("\\n=============================================\\n")
  cat("Session ended:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n")
  cat("Log file:", log_session$file, "\\n")
  cat("=============================================\\n")
  
  sink()
  close(log_session$con)
}

# -----------------------------------------------------------------------------
# EXPORT FUNCTIONS
# -----------------------------------------------------------------------------

# Export table for publication
export_table <- function(data, filename, format = "csv", 
                         digits = 3, na_string = "NA") {
  
  formats <- list(
    csv = function(x, file) write.csv(x, file, row.names = FALSE, na = na_string),
    tsv = function(x, file) write.table(x, file, sep = "\\t", 
                                         row.names = FALSE, na = na_string),
    xlsx = function(x, file) openxlsx::write.xlsx(x, file, 
                                                   asTable = TRUE)
  )
  
  if (!format %in% names(formats)) {
    stop("Unsupported format:", format)
  }
  
  # Round numeric columns
  if (is.data.frame(data)) {
    data <- data %>%
      mutate(across(where(is.numeric), ~ round(., digits)))
  }
  
  # Export
  formats[[format]](data, filename)
  cat("Table exported to:", filename, "\\n")
}

# Create analysis report
create_report <- function(results, template = NULL, output_file = "report.html") {
  if (is.null(template)) {
    # Default report template
    template <- "
# Analysis Report
Generated: {{date}}

## Summary
{{summary}}

## Key Findings
{{findings}}

## Next Steps
{{next_steps}}
"
  }
  
  # Fill template
  report <- template %>%
    str_replace("\\{\\{date\\}\\}", format(Sys.time(), "%Y-%m-%d %H:%M:%S")) %>%
    str_replace("\\{\\{summary\\}\\}", "Analysis completed successfully") %>%
    str_replace("\\{\\{findings\\}\\}", "Check results directory for detailed outputs") %>%
    str_replace("\\{\\{next_steps\\}\\}", "Review figures and statistical tests")
  
  writeLines(report, output_file)
  cat("Report created:", output_file, "\\n")
}
'

writeLines(utility_script, "scripts/utilities/analysis_helpers.R")
cat("   ✓ Created scripts/utilities/analysis_helpers.R\n")

# =============================================================================
# FINAL SUMMARY
# =============================================================================

end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "mins")

cat("\n" + strrep("=", 70) + "\n")
cat("SETUP COMPLETE!\n")
cat(strrep("=", 70) + "\n\n")

cat("Summary:\n")
cat("--------\n")
cat("• R version:", as.character(r_version), "\n")
cat("• CRAN packages:", length(cran_result$installed), "installed\n")
cat("• Bioconductor packages:", length(bioc_result$installed), "installed\n")
cat("• GitHub packages:", length(github_result$installed), "installed\n")
cat("• Project structure created with", length(dirs_to_create), "directories\n")
cat("• Setup duration:", round(duration, 1), "minutes\n\n")

cat("Next steps:\n")
cat("-----------\n")
cat("1. Restart R/RStudio to load the new .Rprofile\n")
cat("2. Run the example script: source('scripts/example_analysis.R')\n")
cat("3. Check the configuration: source('config/settings.R')\n")
cat("4. Review package versions: package_versions.csv\n\n")

cat("Available commands after restart:\n")
cat("---------------------------------\n")
cat("• setup_help()     - Show available helper functions\n")
cat("• check_env()      - Check environment setup\n")
cat("• update_packages()- Update all packages\n")
cat("• clear_temp()     - Clear temporary files\n\n")

cat("Analysis modules:\n")
cat("-----------------\n")
cat("• Morphological analysis: astrocyte_morphological_analysis.Rmd\n")
cat("• Calcium imaging: calcium_imaging_analysis.Rmd\n")
cat("• Cytotoxicity assay: ldh_cytotoxicity_analysis.Rmd\n")
cat("• Transcriptomics: rnaseq_analysis_pipeline.Rmd\n\n")

cat("Happy analyzing! 🧬🔬📊\n")
cat(strrep("=", 70) + "\n")

# Save session info
sessionInfo_file <- paste0("session_info_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
sink(sessionInfo_file)
sessionInfo()
sink()

cat("\nSession info saved to:", sessionInfo_file, "\n")
