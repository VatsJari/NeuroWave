
``` r
---
title: "Transcriptomic Analysis of Astrocytes Post-Microwave Stimulation"
subtitle: "Differential Expression and Pathway Analysis Pipeline"
author: "Your Name/Affiliation"
date: "`r format(Sys.Date(), '%B %d, %Y')`"
output:
  html_document:
    toc: true
    toc_float: true
    toc_depth: 3
    theme: flatly
    highlight: tango
    code_folding: show
    fig_width: 10
    fig_height: 6
    df_print: paged
---
```

```{r setup, include=FALSE, warning=FALSE, message=FALSE}
knitr::opts_chunk$set(
  echo = TRUE,
  warning = FALSE,
  message = FALSE,
  cache = FALSE,
  fig.align = "center"
)
```

# Project Overview

This pipeline performs transcriptomic analysis of astrocytes following
microwave stimulation using RNA-seq data. The analysis includes
differential expression analysis, visualization, and pathway enrichment.

## Analysis Workflow

1.  **Data Preprocessing**: Quality control and normalization
2.  **Differential Expression**: limma-voom analysis
3.  **Visualization**: Volcano plots and expression heatmaps
4.  **Pathway Analysis**: GSEA using MSigDB databases
5.  **Functional Interpretation**: Biological pathway interpretation

# Setup and Installation

## Package Management

```{r package-setup, include=TRUE, warning=FALSE, message=FALSE}
# Load required packages
required_packages <- c(
  # Core analysis
  "tidyverse", "dplyr", "tidyr", "readr",
  # Differential expression
  "limma", "edgeR",
  # Visualization
  "ggplot2", "ggrepel", "viridis", "RColorBrewer", "pheatmap",
  # Pathway analysis
  "fgsea", "msigdbr", "clusterProfiler", "enrichplot",
  # Utilities
  "gridExtra", "cowplot", "ComplexHeatmap", "circlize"
)

# Install missing packages
missing_packages <- required_packages[!(required_packages %in% installed.packages()[ , "Package"])]
if(length(missing_packages)) {
  install.packages(missing_packages)
  
  # Install Bioconductor packages
  if(!"BiocManager" %in% installed.packages()[ , "Package"]) {
    install.packages("BiocManager")
  }
  BiocManager::install(c("limma", "edgeR", "fgsea", "clusterProfiler", 
                        "enrichplot", "ComplexHeatmap"))
}

# Load all packages
invisible(lapply(required_packages, library, character.only = TRUE))

# Set global options
options(stringsAsFactors = FALSE)
set.seed(42)

# Custom color palettes
volcano_colors <- c("Up-regulated" = "#FD841F", 
                   "Down-regulated" = "#38E54D", 
                   "Stable" = "grey90")
gsea_colors <- viridis::viridis(10)
condition_colors <- list(
  "Control" = "#1F77B4",
  "Heat" = "#FF7F0E",
  "Stimulation" = "#E50000",
  "Sham" = "darkgrey"
)
```

# Data Loading and Preprocessing

## Configuration

```{r configuration, include=TRUE}
# Define paths
work_dir <- "/Users/vatsaljariwala/Documents/DIAQNOS/VJ_Microwave_EXP/RNA_seq_thermal_data/RNA"
results_dir <- "~/Desktop/Research Projects (3DBM)/Vatsal/Results/"

# Create results directory if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# Set working directory
setwd(work_dir)
```

## Data Loading Functions

```{r data-loading-functions, include=TRUE}
# Function to load count data
load_count_data <- function(count_file, metadata_file, sample_col = "Sample", 
                           condition_col = "Condition") {
  
  # Load count matrix
  counts <- read.csv(count_file, row.names = 1)
  
  # Load metadata
  metadata <- read.csv(metadata_file)
  
  # Validate data consistency
  if (!all(colnames(counts) %in% metadata[[sample_col]])) {
    warning("Some samples in count matrix are not in metadata")
  }
  
  # Reorder metadata to match count matrix
  metadata <- metadata[match(colnames(counts), metadata[[sample_col]]), ]
  
  # Ensure condition is factor
  metadata[[condition_col]] <- factor(metadata[[condition_col]])
  
  return(list(
    counts = counts,
    metadata = metadata
  ))
}

# Function to perform basic quality control
perform_qc <- function(count_data, min_cpm = 1, min_samples = 2) {
  
  # Create DGEList object
  dge <- DGEList(counts = count_data$counts, 
                 group = count_data$metadata$Condition)
  
  # Calculate library sizes
  dge$samples$lib.size <- colSums(dge$counts)
  
  # Calculate normalization factors
  dge <- calcNormFactors(dge, method = "TMM")
  
  # Filter low-expressed genes
  keep <- rowSums(cpm(dge) >= min_cpm) >= min_samples
  dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]
  
  cat("Original genes:", nrow(dge), "\n")
  cat("Genes after filtering:", nrow(dge_filtered), "\n")
  cat("Percentage retained:", round(nrow(dge_filtered)/nrow(dge)*100, 1), "%\n")
  
  return(dge_filtered)
}

# Function to create design matrix
create_design_matrix <- function(metadata, formula = ~ 0 + Condition) {
  design <- model.matrix(formula, data = metadata)
  colnames(design) <- gsub("Condition", "", colnames(design))
  return(design)
}
```

# Differential Expression Analysis

## Limma-Voom Analysis

```{r limma-analysis-functions, include=TRUE}
# Perform differential expression analysis using limma-voom
perform_limma_analysis <- function(dge_object, design_matrix, contrast_string) {
  
  # Voom transformation
  v <- voom(dge_object, design_matrix, plot = FALSE)
  
  # Fit linear model
  fit <- lmFit(v, design_matrix)
  
  # Define contrasts
  contrasts <- makeContrasts(contrasts = contrast_string, levels = colnames(coef(fit)))
  
  # Fit contrasts
  fit_contrasts <- contrasts.fit(fit, contrasts)
  fit_contrasts <- eBayes(fit_contrasts)
  
  # Extract results
  results <- topTable(fit_contrasts, number = Inf, sort.by = "P")
  
  # Add gene symbols if available
  if (!is.null(rownames(results))) {
    results$SYMBOL <- rownames(results)
    results <- results %>% select(SYMBOL, everything())
  }
  
  return(list(
    results = results,
    voom_object = v,
    fit_object = fit_contrasts,
    contrasts = contrasts
  ))
}

# Filter and annotate results
process_de_results <- function(de_results, pval_threshold = 0.05, 
                              lfc_threshold = 1, adjust_method = "BH") {
  
  # Calculate adjusted p-values if not present
  if (!"adj.P.Val" %in% colnames(de_results)) {
    de_results$adj.P.Val <- p.adjust(de_results$P.Value, method = adjust_method)
  }
  
  # Calculate -log10 p-values
  de_results$neg_log10_pval <- -log10(de_results$P.Value)
  
  # Classify expression changes
  de_results$DE_Status <- case_when(
    de_results$adj.P.Val < pval_threshold & de_results$logFC > lfc_threshold ~ "Up-regulated",
    de_results$adj.P.Val < pval_threshold & de_results$logFC < -lfc_threshold ~ "Down-regulated",
    TRUE ~ "Stable"
  )
  
  # Add significance stars
  de_results$Significance <- case_when(
    de_results$adj.P.Val < 0.001 ~ "***",
    de_results$adj.P.Val < 0.01 ~ "**",
    de_results$adj.P.Val < 0.05 ~ "*",
    TRUE ~ ""
  )
  
  return(de_results)
}
```

## Visualization Functions

```{r visualization-functions, include=TRUE}
# Create volcano plot
create_volcano_plot <- function(de_results, top_n = 10, 
                               pval_threshold = 0.05, lfc_threshold = 1,
                               title = "Differential Expression") {
  
  # Prepare data
  plot_data <- de_results %>%
    mutate(
      highlight = case_when(
        adj.P.Val < pval_threshold & abs(logFC) > lfc_threshold ~ TRUE,
        TRUE ~ FALSE
      )
    )
  
  # Select top genes for labeling
  top_genes <- plot_data %>%
    filter(highlight) %>%
    arrange(adj.P.Val, desc(abs(logFC))) %>%
    head(top_n * 2) %>%  # Get top n up and down
    mutate(
      label_color = ifelse(logFC > 0, volcano_colors[1], volcano_colors[2])
    )
  
  # Create plot
  p <- ggplot(plot_data, aes(x = logFC, y = neg_log10_pval, 
                            color = DE_Status, alpha = highlight)) +
    geom_point(size = 1.5) +
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), 
               linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(pval_threshold), 
               linetype = "dashed", alpha = 0.5) +
    scale_color_manual(values = volcano_colors) +
    scale_alpha_manual(values = c(0.3, 0.8), guide = "none") +
    labs(
      title = title,
      x = "log₂ Fold Change",
      y = "-log₁₀(p-value)",
      color = "Expression"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
  
  # Add gene labels if available
  if (nrow(top_genes) > 0 && "SYMBOL" %in% colnames(top_genes)) {
    p <- p + ggrepel::geom_text_repel(
      data = top_genes,
      aes(label = SYMBOL, color = DE_Status),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.5,
      segment.color = "grey50",
      segment.alpha = 0.5,
      show.legend = FALSE
    )
  }
  
  return(p)
}

# Create enhanced volcano plot with expression levels
create_enhanced_volcano <- function(de_results, top_n = 10) {
  
  p <- ggplot(de_results, aes(x = logFC, y = neg_log10_pval, color = AveExpr)) +
    geom_point(size = 1.5, alpha = 0.7) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
    scale_color_viridis_c(option = "C", name = "Mean\nexpression") +
    labs(
      title = "Differential Expression (colored by expression level)",
      x = "log₂ Fold Change",
      y = "-log₁₀(p-value)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right"
    )
  
  # Add top gene labels
  top_genes <- de_results %>%
    filter(DE_Status != "Stable") %>%
    arrange(adj.P.Val, desc(abs(logFC))) %>%
    head(top_n)
  
  if (nrow(top_genes) > 0 && "SYMBOL" %in% colnames(top_genes)) {
    p <- p + ggrepel::geom_text_repel(
      data = top_genes,
      aes(label = SYMBOL),
      size = 3,
      color = "black",
      box.padding = 0.5,
      max.overlaps = 15
    )
  }
  
  return(p)
}

# Create MA plot
create_ma_plot <- function(de_results) {
  
  p <- ggplot(de_results, aes(x = AveExpr, y = logFC, color = DE_Status)) +
    geom_point(size = 1, alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = "solid", alpha = 0.5) +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
    scale_color_manual(values = volcano_colors) +
    labs(
      title = "MA Plot",
      x = "Average Expression",
      y = "log₂ Fold Change",
      color = "Expression"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right"
    )
  
  return(p)
}
```

# Pathway Analysis

## MSigDB Database Preparation

```{r msigdb-preparation, include=TRUE}
# Prepare MSigDB gene sets
prepare_msigdb_genesets <- function(species = "Homo sapiens") {
  
  # Get available categories
  msigdb_categories <- msigdbr::msigdbr_collections()
  
  # Define which databases to use
  databases <- list(
    Hallmark = list(category = "H", subcategory = NULL),
    GOBP = list(category = "C5", subcategory = "GO:BP"),
    KEGG = list(category = "C2", subcategory = "CP:KEGG"),
    Reactome = list(category = "C2", subcategory = "CP:REACTOME"),
    WikiPathways = list(category = "C2", subcategory = "CP:WIKIPATHWAYS")
  )
  
  # Load each database
  genesets_list <- list()
  
  for (db_name in names(databases)) {
    cat("Loading", db_name, "database...\n")
    
    db_config <- databases[[db_name]]
    geneset <- msigdbr(species = species, 
                      category = db_config$category,
                      subcategory = db_config$subcategory)
    
    # Format for clusterProfiler
    genesets_list[[db_name]] <- geneset %>%
      dplyr::select(gs_name, gene_symbol) %>%
      mutate(database = db_name)
  }
  
  # Combine all databases
  all_genesets <- bind_rows(genesets_list)
  
  # Create list format for fgsea
  geneset_list <- split(all_genesets$gene_symbol, all_genesets$gs_name)
  
  return(list(
    databases = databases,
    genesets = all_genesets,
    geneset_list = geneset_list,
    by_database = split(all_genesets, all_genesets$database)
  ))
}

# Example usage (commented out)
# msigdb_data <- prepare_msigdb_genesets()
```

## GSEA Analysis Functions

```{r gsea-functions, include=TRUE}
# Prepare ranked gene list for GSEA
prepare_ranked_list <- function(de_results, rank_by = "logFC", 
                               gene_col = "SYMBOL") {
  
  # Ensure required columns exist
  required_cols <- c(gene_col, rank_by)
  if (!all(required_cols %in% colnames(de_results))) {
    stop("Missing required columns in DE results")
  }
  
  # Remove NA values
  de_clean <- de_results %>%
    drop_na(all_of(required_cols))
  
  # Create ranked list
  ranked_list <- de_clean[[rank_by]]
  names(ranked_list) <- de_clean[[gene_col]]
  
  # Sort in decreasing order
  ranked_list <- sort(ranked_list, decreasing = TRUE)
  
  return(ranked_list)
}

# Perform GSEA analysis
perform_gsea_analysis <- function(ranked_list, geneset_list, 
                                 min_size = 10, max_size = 500,
                                 nperm = 10000) {
  
  # Run GSEA
  gsea_results <- fgsea::fgsea(
    pathways = geneset_list,
    stats = ranked_list,
    minSize = min_size,
    maxSize = max_size,
    nperm = nperm,
    nproc = 1
  )
  
  # Sort by NES and p-value
  gsea_results <- gsea_results %>%
    arrange(desc(NES), pval)
  
  return(gsea_results)
}

# Perform GSEA for multiple databases
perform_multi_database_gsea <- function(ranked_list, msigdb_data, 
                                       min_size = 10, max_size = 500) {
  
  gsea_results <- list()
  
  for (db_name in names(msigdb_data$by_database)) {
    cat("Running GSEA for", db_name, "...\n")
    
    # Get gene sets for this database
    db_genesets <- msigdb_data$by_database[[db_name]]
    db_list <- split(db_genesets$gene_symbol, db_genesets$gs_name)
    
    # Run GSEA
    tryCatch({
      gsea_res <- perform_gsea_analysis(ranked_list, db_list, 
                                       min_size, max_size)
      gsea_res$database <- db_name
      gsea_results[[db_name]] <- gsea_res
      
      cat("  Found", sum(gsea_res$padj < 0.05), "significant pathways\n")
    }, error = function(e) {
      cat("  Error in", db_name, ":", e$message, "\n")
    })
  }
  
  # Combine all results
  combined_results <- bind_rows(gsea_results)
  
  return(combined_results)
}

# Filter and process GSEA results
process_gsea_results <- function(gsea_results, pval_threshold = 0.05,
                                top_n_per_db = 10) {
  
  # Filter by adjusted p-value
  significant <- gsea_results %>%
    filter(padj < pval_threshold) %>%
    arrange(desc(abs(NES)), padj)
  
  # Get top pathways per database
  top_pathways <- significant %>%
    group_by(database) %>%
    arrange(desc(abs(NES))) %>%
    slice_head(n = top_n_per_db) %>%
    ungroup()
  
  # Add direction classification
  top_pathways <- top_pathways %>%
    mutate(
      direction = ifelse(NES > 0, "Up-regulated", "Down-regulated"),
      pathway_name = gsub("_", " ", pathway),
      pathway_name = str_to_title(pathway_name)
    )
  
  return(list(
    all_results = gsea_results,
    significant = significant,
    top_pathways = top_pathways
  ))
}
```

## GSEA Visualization Functions

```{r gsea-visualization, include=TRUE}
# Create GSEA dot plot
create_gsea_dotplot <- function(gsea_results, top_n = 15, 
                               title = "GSEA Results") {
  
  # Prepare data
  plot_data <- gsea_results %>%
    filter(padj < 0.05) %>%
    arrange(desc(abs(NES))) %>%
    head(top_n) %>%
    mutate(
      pathway_label = gsub("_", " ", pathway),
      pathway_label = str_wrap(pathway_label, width = 40),
      pathway_label = factor(pathway_label, 
                            levels = rev(pathway_label))
    )
  
  p <- ggplot(plot_data, aes(x = NES, y = pathway_label, 
                            size = -log10(padj), color = NES)) +
    geom_point(alpha = 0.8) +
    scale_size_continuous(range = c(2, 6), 
                         name = "-log₁₀(adj.p)") +
    scale_color_gradient2(
      low = "#2E86AB",
      mid = "white",
      high = "#A23B72",
      midpoint = 0,
      name = "NES"
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    labs(
      title = title,
      x = "Normalized Enrichment Score (NES)",
      y = "Pathway"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 9),
      panel.grid.major.y = element_blank(),
      legend.position = "right"
    )
  
  return(p)
}

# Create GSEA enrichment plot for specific pathway
create_enrichment_plot <- function(ranked_list, geneset_list, 
                                  pathway_name, title = NULL) {
  
  if (!pathway_name %in% names(geneset_list)) {
    stop("Pathway not found in geneset list")
  }
  
  # Create plot
  p <- enrichplot::gseaplot2(
    geneList = ranked_list,
    geneSet = pathway_name,
    pvalue_table = TRUE,
    title = ifelse(is.null(title), pathway_name, title)
  )
  
  return(p)
}

# Create multi-database GSEA summary plot
create_gsea_summary_plot <- function(gsea_processed, top_n_per_db = 5) {
  
  # Get top pathways from each database
  top_by_db <- gsea_processed$top_pathways %>%
    group_by(database) %>%
    slice_head(n = top_n_per_db) %>%
    ungroup()
  
  # Create faceted plot
  p <- ggplot(top_by_db, aes(x = NES, y = reorder(pathway, NES), 
                            fill = database, size = -log10(padj))) +
    geom_point(shape = 21, alpha = 0.8) +
    facet_wrap(~ database, scales = "free_y", ncol = 2) +
    scale_fill_brewer(palette = "Set2", name = "Database") +
    scale_size_continuous(range = c(1, 4), name = "-log₁₀(adj.p)") +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    labs(
      title = "Top Enriched Pathways by Database",
      x = "Normalized Enrichment Score (NES)",
      y = "Pathway"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 8),
      panel.grid.major.y = element_blank()
    )
  
  return(p)
}
```

# Complete Analysis Pipeline

## Main Analysis Function

```{r main-analysis-function, include=TRUE}
run_rna_seq_analysis <- function(
  count_file,
  metadata_file,
  contrast = "Heat - Control",
  output_dir = "./results",
  pval_threshold = 0.05,
  lfc_threshold = 1,
  species = "Homo sapiens"
) {
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat("Starting RNA-seq analysis pipeline...\n")
  cat("Output directory:", output_dir, "\n")
  
  # Step 1: Load data
  cat("\n1. Loading data...\n")
  data <- load_count_data(count_file, metadata_file)
  cat("   Samples:", ncol(data$counts), "\n")
  cat("   Genes:", nrow(data$counts), "\n")
  cat("   Conditions:", toString(levels(data$metadata$Condition)), "\n")
  
  # Step 2: Quality control
  cat("\n2. Performing quality control...\n")
  dge_qc <- perform_qc(data)
  
  # Step 3: Create design matrix
  cat("\n3. Creating design matrix...\n")
  design <- create_design_matrix(data$metadata)
  cat("   Design matrix dimensions:", dim(design), "\n")
  
  # Step 4: Differential expression analysis
  cat("\n4. Performing differential expression analysis...\n")
  de_results <- perform_limma_analysis(dge_qc, design, contrast)
  
  # Process DE results
  de_processed <- process_de_results(
    de_results$results,
    pval_threshold = pval_threshold,
    lfc_threshold = lfc_threshold
  )
  
  # Summary statistics
  de_summary <- de_processed %>%
    count(DE_Status) %>%
    mutate(Percentage = round(n / nrow(de_processed) * 100, 1))
  
  cat("   DE summary:\n")
  print(de_summary)
  
  # Step 5: Prepare MSigDB databases
  cat("\n5. Preparing MSigDB databases...\n")
  msigdb_data <- prepare_msigdb_genesets(species = species)
  
  # Step 6: Prepare ranked list for GSEA
  cat("\n6. Preparing ranked gene list for GSEA...\n")
  ranked_list <- prepare_ranked_list(de_processed)
  
  # Step 7: Perform GSEA
  cat("\n7. Performing GSEA analysis...\n")
  gsea_results <- perform_multi_database_gsea(ranked_list, msigdb_data)
  
  # Process GSEA results
  gsea_processed <- process_gsea_results(gsea_results)
  
  # GSEA summary
  gsea_summary <- gsea_processed$significant %>%
    group_by(database) %>%
    summarise(
      n_sig = n(),
      max_NES = round(max(abs(NES)), 2),
      .groups = "drop"
    )
  
  cat("   GSEA summary:\n")
  print(gsea_summary)
  
  # Step 8: Create visualizations
  cat("\n8. Creating visualizations...\n")
  
  # Volcano plots
  volcano1 <- create_volcano_plot(
    de_processed,
    title = paste("Differential Expression:", contrast)
  )
  
  volcano2 <- create_enhanced_volcano(de_processed)
  
  # MA plot
  ma_plot <- create_ma_plot(de_processed)
  
  # GSEA plots
  gsea_dotplot <- create_gsea_dotplot(gsea_results)
  gsea_summary_plot <- create_gsea_summary_plot(gsea_processed)
  
  # Step 9: Save results
  cat("\n9. Saving results...\n")
  
  # Save DE results
  de_file <- file.path(output_dir, "differential_expression_results.csv")
  write_csv(de_processed, de_file)
  cat("   DE results saved to:", de_file, "\n")
  
  # Save GSEA results
  gsea_file <- file.path(output_dir, "gsea_results.csv")
  write_csv(gsea_results, gsea_file)
  cat("   GSEA results saved to:", gsea_file, "\n")
  
  # Save significant DE genes
  sig_de <- de_processed %>%
    filter(DE_Status != "Stable")
  
  sig_file <- file.path(output_dir, "significant_genes.csv")
  write_csv(sig_de, sig_file)
  cat("   Significant genes saved to:", sig_file, "\n")
  
  # Save plots
  plots <- list(
    volcano_standard = volcano1,
    volcano_enhanced = volcano2,
    ma_plot = ma_plot,
    gsea_dotplot = gsea_dotplot,
    gsea_summary = gsea_summary_plot
  )
  
  for (plot_name in names(plots)) {
    plot_file <- file.path(output_dir, paste0(plot_name, ".pdf"))
    ggsave(plot_file, plots[[plot_name]], width = 8, height = 6)
    cat("   Plot saved:", plot_file, "\n")
  }
  
  # Create summary report
  create_summary_report(de_processed, gsea_processed, output_dir)
  
  cat("\nAnalysis complete!\n")
  
  return(list(
    data = data,
    dge_object = dge_qc,
    de_results = de_processed,
    gsea_results = gsea_processed,
    plots = plots,
    output_dir = output_dir
  ))
}

# Create summary report
create_summary_report <- function(de_results, gsea_results, output_dir) {
  
  report_text <- c(
    "# RNA-seq Analysis Summary",
    paste("Generated on:", Sys.Date()),
    "",
    "## Differential Expression",
    "",
    paste("Total genes analyzed:", nrow(de_results)),
    paste("Up-regulated genes:", sum(de_results$DE_Status == "Up-regulated")),
    paste("Down-regulated genes:", sum(de_results$DE_Status == "Down-regulated")),
    paste("Significant genes (FDR < 0.05):", sum(de_results$adj.P.Val < 0.05)),
    "",
    "## GSEA Results",
    "",
    paste("Total pathways tested:", nrow(gsea_results$all_results)),
    paste("Significant pathways (FDR < 0.05):", nrow(gsea_results$significant)),
    "",
    "### By Database",
    ""
  )
  
  # Add database-specific summaries
  for (db_name in unique(gsea_results$all_results$database)) {
    db_sig <- gsea_results$significant %>%
      filter(database == db_name) %>%
      nrow()
    
    report_text <- c(report_text,
                     paste("-", db_name, ":", db_sig, "significant pathways"))
  }
  
  # Add top pathways
  report_text <- c(report_text,
                   "",
                   "### Top Enriched Pathways",
                   "")
  
  top_pathways <- gsea_results$top_pathways %>%
    arrange(desc(abs(NES))) %>%
    head(10)
  
  for (i in 1:nrow(top_pathways)) {
    pathway <- top_pathways[i, ]
    report_text <- c(report_text,
                     paste(i, ".", pathway$pathway, 
                           "(NES =", round(pathway$NES, 2),
                           ", FDR =", format(pathway$padj, scientific = TRUE, digits = 2), ")"))
  }
  
  # Write report
  report_file <- file.path(output_dir, "analysis_summary.md")
  writeLines(report_text, report_file)
  cat("Summary report saved to:", report_file, "\n")
}
```

## Example Usage

```{r example-usage, eval=FALSE}
# Example usage of the complete pipeline
results <- run_rna_seq_analysis(
  count_file = "counts_astrocytes.csv",
  metadata_file = "samp_desc_astro.csv",
  contrast = "Heat - Control",
  output_dir = "./analysis_results",
  pval_threshold = 0.05,
  lfc_threshold = 1
)

# Access results
print(results$de_results %>% filter(DE_Status != "Stable") %>% head(10))
print(results$gsea_results$top_pathways)

# Display plots
print(results$plots$volcano_standard)
print(results$plots$gsea_dotplot)
```

# Quick Analysis Example

## Direct Analysis Code

```{r direct-analysis, include=TRUE}
# This section provides the direct analysis code from the original script
# in a cleaned and organized format

analyze_astrocytes_direct <- function() {
  
  # Load data
  new_data <- read.csv("counts_astrocytes.csv", row.names = 1)
  metadata <- read.csv("samp_desc_astro.csv")
  
  # Create DGEList
  d0 <- DGEList(new_data)
  d0 <- calcNormFactors(d0, method = "TMM")
  
  # Filter low-expressed genes
  cutoff <- 1
  keep <- rowSums(cpm(d0) >= cutoff) >= 2
  d <- d0[keep, ]
  cat("Genes after filtering:", dim(d), "\n")
  
  # Create design matrix
  design <- model.matrix(~ 0 + Condition, data = metadata)
  colnames(design) <- gsub("Condition", "", colnames(design))
  
  # Log normalize
  log2cpm <- cpm(d0, log = TRUE, prior.count = 1)
  
  # Fit model and compute contrasts
  fit <- lmFit(log2cpm, design)
  contr <- makeContrasts(Heat - Control, levels = colnames(coef(fit)))
  
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  
  # Save results
  write.csv(top.table, file.path(results_dir, "differential_expression.csv"))
  
  # Return results
  return(list(
    de_results = top.table,
    design = design,
    contrasts = contr
  ))
}

# Run the analysis
direct_results <- analyze_astrocytes_direct()
```

# Session Information

```{r session-info, include=TRUE}
sessionInfo()
```

# Usage Instructions

## Quick Start

1.  **Prepare Data**:
    -   Count matrix (CSV with genes as rows, samples as columns)
    -   Metadata file (CSV with sample information)
2.  **Configure Analysis**:
    -   Set file paths in the configuration section
    -   Define contrast of interest (e.g., "Treatment - Control")
3.  **Run Analysis**:
    -   Use `run_rna_seq_analysis()` for comprehensive analysis
    -   Or use individual functions for specific steps
4.  **Review Results**:
    -   Check output directory for results and plots
    -   Review summary report for key findings

## Customization Options

-   **Statistical Thresholds**: Adjust p-value and fold-change
    thresholds
-   **Normalization Methods**: Change TMM to other methods as needed
-   **Pathway Databases**: Select specific MSigDB databases
-   **Visualization Styles**: Modify plot aesthetics in visualization
    functions

## Output Files

The pipeline generates: - Differential expression results - GSEA pathway
enrichment results - Multiple plot formats (PDF) - Summary statistics
report - Processed data files

## Troubleshooting

-   **Memory Issues**: Filter genes more stringently or analyze subsets
-   **Convergence Problems**: Increase permutation count in GSEA
-   **Missing Gene Symbols**: Ensure proper annotation of count matrix
-   **Design Matrix Errors**: Check metadata consistency with count
    matrix

------------------------------------------------------------------------

*This transcriptomic analysis pipeline provides comprehensive
differential expression and pathway analysis for RNA-seq data following
microwave stimulation experiments. For technical support or custom
modifications, please contact the authors.*

**Citation**: If using this pipeline, please cite: Limma (Ritchie et
al., 2015), edgeR (Robinson et al., 2010), and MSigDB (Liberzon et al.,
2015) \`\`\`
