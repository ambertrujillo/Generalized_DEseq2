################################################
##### Preprocessing data for DEGs analysis #####
################################################
work_dir='/path/to/working/directory/'
data='/path/to/data.h5ad'
metadata='/path/to/metadata.csv'

setwd(work_dir)


library(SummarizedExperiment)
library(zellkonverter)
library(harmony)
library(Seurat)
library(SeuratDisk)
library(reticulate)
library(anndata)
library(tidyverse)
library(ggplot2)

##### All clusters #####
# Gather raw pseudobulk data
data <- read_h5ad(data)

# Get count data from .h5ad file
counts = data$layers["counts"]
counts = CreateSeuratObject(counts = t(as.matrix(counts)), meta.data=data$obs)

# Calculate pseudobulk data by SampleID
pseudobulk = AggregateExpression(counts, group.by = "SampleID", assay = "RNA", slot = "counts")

# Extract pseudobulked count data from Seurat object and make it a matrix
count_matrix = pseudobulk$RNA
count_matrix = as.matrix(count_matrix)

# Work with Metadata
### Read in the metdata to see the sample numbers
metadata = read.csv(metadata)

length(unique(metadata$SampleID))


### Fix the naming issues that SeuratAggregate introduced
colnames(count_matrix) <- gsub("^g", "", colnames(count_matrix))
count_matrix = as.matrix(count_matrix)

table(metadata$SampleID %in% colnames(count_matrix))

# See missing samples Count matrix
missing_in_metadata <- setdiff(colnames(count_matrix), metadata$SampleID)
print(missing_in_metadata) # None missing between count matrix and metadata

counts <- count_matrix[, colnames(count_matrix) %in% metadata$SampleID]

missing_in_sc <- setdiff(metadata$SampleID, colnames(counts))
print(missing_in_sc)
# We do not have sequencing data for ("D19-12386" "D19-12393" "Control4" "Control10" "FTLD-GRN1")
# Remove those we do not have sequencing data for 
metadata = subset(metadata, !(SampleID %in% missing_in_sc))

table(metadata$general_disease)

write.csv(metadata, file = "outs/metadata.csv", row.names = FALSE)
write.csv(counts, file = "outs/counts.csv")







