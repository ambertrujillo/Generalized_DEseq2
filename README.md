# Generalized_DESeq2
**Generalized_DESeq2** is an R-based workflow for identifying differentially expressed genes (DEGs) from single-cell RNA sequencing (scRNA-seq) data using the DESeq2 framework.

The pipeline follows the general DESeq2 best-practices workflow while incorporating customized preprocessing steps tailored for single-cell data. Initial quality control and generation of read count matrices are performed in Python prior to running this pipeline.

## Required libraries
- `SummarizedExperiment`
- `zellkonverter`
- `harmony`
- `Seurat`
- `SeuratDisk`
- `reticulate`
- `anndata`
- `tidyverse`
- `ggplot2`
- `DESeq2`
- `airway`
- `magrittr`
- `EnhancedVolcano`

It is recommended to install all dependencies before running the pipeline. Ensure compatible package versions, particularly for Seurat and zellkonverter when working with `.h5ad` files.

## Required inputs
The pipeline requires the following input files:
- Read count matrix
    - Format: `.h5ad`
    - Generated from prior Python-based preprocessing
    - Must contain raw count data suitable for DESeq2 input

- Metadata file
    - Format: `.csv`
    - Must include:
        - Sample identifiers
        - Experimental condition(s)
        - Variable(s) of interest for differential expression analysis

***Note: Ensure that sample identifiers in the metadata file match those in the count matrix.***

## Running the pipeline
The pipeline can be run interactively in R in the following order:
1. Preprocess.R
    - Prepares read count matrix and metadata for differential gene expression analysis.
    - Outputs:
       - A processed count matrix (`outs/counts.csv`)
       - A matched metadata file (`outs/metadata.csv`)
     These files are aligned and formatted for downstream DEG analysis.
2. DEGs.R
    - Runs differential gene expression analysis and visualizes data (MA and volcano plots)
    - This script is designed to filter out genes that have less than 10 counts across at least 3 individuals. To change this, edit lines 59 and 60.
    - There is an option to run with contrasts (line 52-55, 71-73, and 86-87)
    - Outputs:
       - Controls vs. Conditon DEGs before LFCShrink (`outs/preShrink.csv`)
       - Controls vs. Conditon DEGs after LFCShrink (`outs/DEGs.csv`)

The necessary scripts can be found in scripts/
***Note: Be sure to change the paths in the scripts to fit your data***