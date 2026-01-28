# NeuroWave: Exploring Microwave exposure effects on glial cells 

A comprehensive computational pipeline for analyzing cellular responses to microwave stimulation across multiple modalities: morphology, calcium dynamics, cytotoxicity, and transcriptomics.

## 📊 Overview

This repository contains a complete analysis suite for studying cellular responses to microwave stimulation. The pipeline integrates data from:
- **Morphological analysis** (CellProfiler outputs)
- **Calcium imaging** (GCaMP fluorescence traces)
- **Cytotoxicity assays** (LDH measurements)
- **Transcriptomics** (RNA-seq data)

Each analysis module is designed to be modular, reproducible, and produces publication-ready visualizations.

## 📁 Repository Structure

```
microwave_analysis/
├── README.md                     # This file
├── LICENSE                       # License information
├── .gitignore                    # Git ignore file
├── setup_environment.R           # Environment setup script
├── morphological_analysis/       # Cell morphology analysis
│   ├── Astrocyte_morphological_analysis.Rmd
│   ├── Microglia_morphological_analysis.Rmd
│   ├── GBM_morphological_analysis.Rmd
├── calcium_imaging/              # Calcium dynamics analysis
│   ├── Calcium_Signalling_Analysis.Rmd
├── cytotoxicity_assay/           # LDH assay analysis
│   ├── Ldh_Cytotoxicity_Analysis.Rmd
├── transcriptomics/              # RNA-seq analysis
│   ├── Transcriptomics_BulkSeq_Analysis.Rmd
│   ├── msigdb_databases/         # Pathway databases
└──── 

```

## 🚀 Quick Start

### Prerequisites

1. **R (version 4.1.0 or higher)**
2. **RStudio (recommended)**
3. **Required system libraries**:
   - macOS: `brew install libpng cairo`
   - Ubuntu: `sudo apt-get install libcurl4-openssl-dev libssl-dev libxml2-dev`
   - Windows: Install Rtools

### Installation

1. **Clone the repository**:
   ```bash
   git clone https://github.com/yourusername/microwave-analysis.git
   cd microwave-analysis
   ```

2. **Set up R environment**:
   Open R/RStudio and run:
   ```r
   # Run setup script
   source("setup_environment.R")
   
   # Or manually install packages
   install.packages(c("tidyverse", "ggplot2", "ggpubr", "patchwork", "rmarkdown"))
   ```

3. **Test installation**:
   ```r
   # Check if all required packages are installed
   source("utilities/check_installation.R")
   ```

## 📋 Step-by-Step Usage

### 1. Morphological Analysis

**For analyzing cell morphology from CellProfiler outputs:**

#### Step 1: Prepare your data
- Place CellProfiler CSV outputs in a structured folder:
  ```
  your_data/
  ├── primary_objects.csv    # Soma measurements
  └── secondary_objects.csv  # Whole cell measurements
  ```

### 2. Calcium Imaging Analysis

**For analyzing GCaMP calcium traces:**

#### Step 1: Organize calcium trace files
```
calcium_data/
├── experiment1.csv
├── experiment2.csv
└── metadata.csv    # Optional: experimental conditions
```

#### Step 2: Configure analysis parameters
Edit the calcium analysis script:
```r
# Set your parameters
data_folder <- "path/to/your/calcium/data"
results_folder <- "path/to/results"
recording_time_min <- 3        # Recording duration in minutes
cutoff_fixed <- 0.2            # Activity threshold
delta_multiplier <- 2          # Peak detection sensitivity
```


### 3. LDH Cytotoxicity Assay

**For analyzing LDH release data:**

#### Step 1: Prepare your LDH data
Create a CSV file with columns:
- `Condition` (e.g., "G0", "G1", etc.)
- `Lum_score` (raw luminescence values)
- `Replicate` (optional)

#### Step 2: Run the analysis
```r
# Load the analysis script
rmarkdown::render("cytotoxicity_assay/ldh_cytotoxicity_analysis.Rmd",
                  params = list(
                    data_path = "path/to/your/ldh_data.csv",
                    cell_type = "GBM",  # or "Astrocyte" or "Microglia"
                    output_dir = "ldh_results"
                  ))
```

#### Step 3: Interpret results
Output includes:
- Normalized cytotoxicity values
- Statistical comparisons
- Dose-response visualizations
- Comparative analysis across cell types

### 4. Transcriptomics Analysis

**For analyzing RNA-seq data:**

#### Step 1: Prepare count data
- Count matrix (genes × samples)
- Metadata file (samples × conditions)


#### Step 2: Explore results
The pipeline produces:
- Differential expression results
- Volcano plots
- Pathway enrichment analysis
- Interactive visualizations

## 🔧 Configuration Guide

### Experimental Design Mapping

Map your experimental conditions to the analysis codes:

| Condition Code | MW Status | Exposure Time | Recovery Time |
|----------------|-----------|---------------|---------------|
| G0, A0, M0     | ON        | 5 min         | 0 h           |
| G1, A1, M1     | OFF       | 5 min         | 0 h           |
| G2, A2, M2     | ON        | 30 min        | 0 h           |
| G3, A3, M3     | OFF       | 30 min        | 0 h           |
| G4, A4, M4     | ON        | 5 min         | 24 h          |
| G5, A5, M5     | OFF       | 5 min         | 24 h          |
| G6, A6, M6     | ON        | 30 min        | 24 h          |
| G7, A7, M7     | OFF       | 30 min        | 24 h          |

### Customizing Analysis Parameters

Each analysis module has configurable parameters:

#### Morphological Analysis
```r
# In the R Markdown files, adjust:
company_colors <- c("#E50000", "#008A8A", ...)  # Color palette
window_size <- 50                               # Rolling window size
optimal_k <- 6                                  # Number of clusters
pval_threshold <- 0.05                          # Statistical threshold
```

#### Calcium Imaging
```r
# In calcium_imaging_analysis.Rmd:
recording_time_min <- 3        # Recording duration
cutoff_fixed <- 0.2            # Activity threshold
delta_multiplier <- 2          # Peak detection sensitivity
skip_correction <- TRUE        # Skip baseline correction
```

#### Statistical Methods
- **Morphology**: KS tests, Wilcoxon tests, Chi-square tests
- **Calcium**: Peak detection algorithms, firing rate calculations
- **Cytotoxicity**: t-tests, ANOVA, dose-response modeling
- **Transcriptomics**: limma-voom, DESeq2, GSEA

## 📈 Output Interpretation

### Morphological Analysis Outputs

1. **Cell Density Plots**: Normalized cell counts per condition
2. **Parameter Distributions**: Box plots of morphological features
3. **Cluster Analysis**: UMAP/t-SNE visualizations of phenotypes
4. **Phenotype Proportions**: Bar plots of morphological classifications

### Calcium Imaging Outputs

1. **Trace Visualizations**: Raw and processed calcium traces
2. **Event Detection**: Raster plots of calcium events
3. **Firing Rate Analysis**: Statistical comparisons of activity
4. **Waveform Analysis**: Average calcium transient shapes

### Statistical Significance Codes

- \* p < 0.05
- \*\* p < 0.01
- \*\*\* p < 0.001
- ns not significant

## 🔍 Troubleshooting

### Common Issues and Solutions

#### Issue 1: Package Installation Failures
```r
# If Bioconductor packages fail to install:
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.16")  # Specify version
```

#### Issue 2: Memory Errors
```r
# Increase memory limit
memory.limit(size = 16000)  # On Windows

# Use 64-bit R
# Filter data early in pipeline
# Process in chunks for large datasets
```

#### Issue 3: File Path Errors
```r
# Use relative paths
data_path <- file.path("data", "your_file.csv")

# Check file existence
if (!file.exists(data_path)) {
  stop("File not found: ", data_path)
}

# Set working directory properly
setwd("/full/path/to/your/project")
```

#### Issue 4: Missing Dependencies
```r
# Check and install missing packages
required <- c("tidyverse", "ggplot2", "dplyr")
missing <- required[!(required %in% installed.packages()[,"Package"])]
if(length(missing)) install.packages(missing)
```

### Debugging Tips

1. **Run in stages**: Test each analysis module separately
2. **Check data formats**: Ensure CSV files are properly formatted
3. **Review log files**: Check R console output for warnings/errors
4. **Use example data**: Test with provided example datasets first


## 🤝 Contributing

We welcome contributions! Please follow these steps:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

### Code Style Guidelines

- Use meaningful variable names
- Comment complex sections
- Follow tidyverse style guide
- Include example data for new features
- Update documentation accordingly

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📞 Support

For questions, issues, or feature requests:

1. **Check the documentation**: Review this README and inline documentation
2. **Search existing issues**: Your question might already be answered
3. **Open a new issue**: Provide a minimal reproducible example
4. **Contact maintainers**: [Vatsal D. Jariwala]

## 📊 Citation

If you use this software in your research, please cite:

```bibtex
@software{microwave_analysis_suite,
  author = Vatsal D. Jariwala},
  title = {Microwave Stimulation Analysis Suite},
  year = {2025},
  url = {https://github.com/VatsJari/NeuroWave}
}
```

## 🎯 Key Features

- **Modular Design**: Each analysis module works independently
- **Reproducible**: Complete documentation and version control
- **Scalable**: Handles datasets from small to large scale
- **Customizable**: Easily adjust parameters for specific needs
- **Visualization-Rich**: Publication-ready plots and figures
- **Statistical Rigor**: Multiple testing corrections and effect sizes

---

**Happy Analyzing!** 🧬🔬📊

*For best results, start with the example datasets and gradually replace them with your own data. Don't hesitate to open issues if you encounter problems or have suggestions for improvements.*
