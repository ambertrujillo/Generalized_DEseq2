###############################################################################
##### Run DEGs on ExNeu1 for disease with Control as the Reference #####
###############################################################################
work_dir='/path/to/working/directory/'

setwd(work_dir)

library(DESeq2)
library(ggplot2)
library(airway)
library(magrittr)
library(EnhancedVolcano)

##### Make sure that your count matrix matches the order of your metadata
metadata = read.csv("outs/metadata.csv")

## Remove unnecessary column introduced when saving file in Preprocess.R
counts = read.csv("outs/counts.csv")
rownames(counts) = counts$X
counts = subset(counts, select=-X)

## Fix the sample names so the match and are compatible with DEseq2
colnames(counts) = gsub("\\.", "_", colnames(counts))
colnames(counts) = gsub("X", "", colnames(counts))
metadata$SampleID = gsub("-", "_", metadata$SampleID)

### See if SampleIDs are in colnames and if they match
table(colnames(counts) %in% metadata$SampleID)

# Change the rownames of the metadata to sample names
rownames(metadata) = metadata$SampleID #Change the rownames of the meta data to the sample names
all(rownames(metadata) == colnames(counts))
counts <- counts[, rownames(metadata)]
all(rownames(metadata) == colnames(counts))

# Turning general_disease and study into a factor
metadata$general_disease = as.factor(metadata$general_disease)
metadata$study = as.factor(metadata$study)


### Running DEseq2
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ <covariates>)

# Replace `<covariates>` with the variable(s) of interest (e.g., `condition`, `diagnosis`, or additional covariates such as `age + sex + condition`) according to your experimental design.

# Relevel so Control is the reference
dds$<covariate-of-interest> <- relevel(dds$general_disease, ref = "Control")

## TO RUN WITH CONTRASTS (No need to relevel since you are specifying the contrasts)
# dds <- DESeqDataSetFromMatrix(countData = counts,
#                              colData = metadata,
#                              design = ~ 0 + <covariates>)


# filter out genes that have less than 10 counts across at least 3 individuals.
individual_counts <- rowSums(counts(dds) >= 10)
keep <- individual_counts >= 3
dds <- dds[keep,]

### Run DEseq2
dds = DESeq(dds)
# Get results names for contrasts
resultsNames(dds)

noShrink_res <- results(dds, independentFiltering = FALSE, name = "name_of_comparison_of_interest_listed_in_resultsNames(dds)")
noShrink_res

## TO RUN WITH CONTRASTS - In "contrast = " the first position = covariate of interest, second positon = treatment condition, third position = control
#noShrink_res <- results(dds, independentFiltering = FALSE, contrast = c("covariate", "treatment", "control"))
#noShrink_res


write.csv(as.data.frame(noShrink_res), 
          file="outs/preShrink.csv")

# Shrinkage for visualization and ranking
# Get names for specification of results desired
resultsNames(dds)

# Shrinnk results for both comparisons: FTD vs. Control and AD vs. Control
shrink_res = lfcShrink(dds, coef="name_of_comparison_of_interest_listed_in_resultsNames", type="apeglm", res = noShrink_res)

## TO RUN WITH CONTRASTS - In "contrast = " the first position = covariate of interest, second positon = treatment condition, third position = control
#shrink_res = lfcShrink(dds, contrast = c("covariate", "treatment", "control"), type="ashr", res = noShrink_res)


# Reorder results by smalles p-value
resOrdered <- shrink_res[order(shrink_res$pvalue),]

# How many adj p-valuyes are less than 0.1
sum(resOrdered$padj < 0.05, na.rm=TRUE)


write.csv(as.data.frame(resOrdered), 
          file="outs/DEGs.csv")


##### Visualize #####
plotMA(noShrink, ylim=c(-2,2))
plotMA(resOrdered, ylim=c(-2,2))


### Volcano plot #####
# Create dataframe for plotting
res_df <- data.frame(resOrdered)
res_df$gene <- rownames(res_df)

##### Enhanced volcano #####
library(airway)
library(magrittr)
library(EnhancedVolcano)

EnhancedVolcano(res_df,
                lab = rownames(res_df),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 6.0,
                drawConnectors = TRUE,
                ylab = '-Log10(padj)',
                subtitle = "Differential Gene Expression")

