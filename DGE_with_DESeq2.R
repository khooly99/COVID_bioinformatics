# Install required packages (if not already installed)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  BiocManager::install("tidyverse")
}
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}
if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
  BiocManager::install("AnnotationDbi")
}
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  BiocManager::install("pheatmap")
}
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  BiocManager::install("RColorBrewer")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

# Load packages
library(DESeq2)
library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)


## -----------------------------------------------------------------------------
# Import raw gene counts
count <- read.table("data/count.txt", header=TRUE, row.names=1)
head(count)

# Import metadata
meta <- read.csv(file="data/meta.csv", header=TRUE)
head(meta)


## -----------------------------------------------------------------------------
# Pre-processing of counts and metadata
# Rotate metadata
meta <- t(meta)

# Check classes of count and metadata
class(count)
class(meta)

# (If required) Change to dataframe 
meta <- as.data.frame(meta)
count <- as.data.frame(count)

# Check if number of samples are equal in count and metadata
table(colnames(count) == rownames(meta))

# Standardise values in columns of positive or negative covid
meta$Sample_characteristics_ch[meta$Sample_characteristics_ch == "sars-cov-2 positivity: pos"] <- "positive"
meta$Sample_characteristics_ch[meta$Sample_characteristics_ch == "sars-cov-2 positivity: neg"] <- "negative"

# (If required) Select only COVID and controls
meta <- subset(meta, Sample_characteristics_ch %in% c("positive", "negative"))
match <- match(rownames(meta), colnames(count))
count <- count[, match]
table(colnames(count) == rownames(meta))

## -----------------------------------------------------------------------------
# Creating DESeq2 Dataset
# (If required) Add batch as design
dds <- DESeqDataSetFromMatrix(countData = count, 
                              colData=meta, 
                              design= ~Sample_characteristics_batch + Sample_characteristics_ch)
dim(dds)

# Filtering low count genes
summary(dds$Sample_characteristics_ch == "positive")
summary(dds$Sample_characteristics_ch == "negative")
# Obtain lowest number of samples from positive or negative (in this case = 93)
keep <- rowSums(counts(dds) >= 10) >= 93
dds <- dds[keep,]

# Set up factor levels
dds$Sample_characteristics_ch <- relevel(dds$Sample_characteristics_ch, ref = "negative")
dds$Sample_characteristics_ch
dds$Sample_characteristics_batch <- factor(dds$Sample_characteristics_batch)
dds$Sample_characteristics_batch

## -----------------------------------------------------------------------------
# DESeq2 Analysis
dds <- DESeq(dds, parallel=TRUE)

## -----------------------------------------------------------------------------
# Transformation
vsd <- vst(dds, blind=FALSE)
vsd_mat <- assay(vsd)

## -----------------------------------------------------------------------------
# Plots for QC
# PCA plot
plotPCA(vsd, intgroup=c("Sample_characteristics_ch"))
plotPCA(vsd, intgroup=c("Sample_characteristics_batch"))

# Batch-corrected PCA
batch_corrected <- limma::removeBatchEffect(assay(vsd),
                                            batch = dds$Sample_characteristics_batch)
# Visualise batch-corrected PCA
assay(vsd) <- batch_corrected
plotPCA(vsd, intgroup=c("Sample_characteristics_ch"))
plotPCA(vsd, intgroup=c("Sample_characteristics_batch"))

# Save vsd values
write.table(batch_corrected, "vsd.txt")

## -----------------------------------------------------------------------------
# DESeq2 result
res05 <- results(dds, alpha = 0.05)
res05 <- na.omit(res05)

# Order result based on p-adj
res05_DF <- res05[order(res05$padj),]
head(as.data.frame(res05_DF))

# Find significant genes 
res_sig_up <- res05_DF[which(res05_DF$padj < 0.05 & res05_DF$log2FoldChange > 2),]
dim(res_sig_up)
res_sig_down <- res05_DF[which(res05_DF$padj < 0.05 & res05_DF$log2FoldChange < -2),]
dim(res_sig_down)

# Add columns for labels of DEGs (NO: non-DEGs, UP: upregulated, DOWN: downregulated genes)
res05_DF$DGE <- "NO"
res05_DF$DGE[res05_DF$log2FoldChange > 2 & res05_DF$padj < 0.05] <- "UP"
res05_DF$DGE[res05_DF$log2FoldChange < -2 & res05_DF$padj < 0.05] <- "DOWN"

## -----------------------------------------------------------------------------
# (If required) Annotation of results from ENSEMBL to GENE SYMBOLS
annotations_orgDb <- AnnotationDbi::select(org.Hs.eg.db, # database
                                           keys = res05_DF$GeneName,  # data to use for retrieval
                                           columns = "SYMBOL", # information to retrieve for given data
                                           keytype = "ENSEMBL") # type of data given in 'keys' argument
colnames(res05_DF)[which(names(res05_DF) == "GeneName")] <- "ENSEMBL"
res05_DF <- merge(res05_DF, annotations_orgDb, by.x = , by.y = "ENSEMBL")
res05_DF <- na.omit(res05_DF)
summary(res05_DF)

# Export results to txt
write.table(res05_DF, file="res05.txt")

## -----------------------------------------------------------------------------
# Visualisation of results
# MA Plot
DESeq2::plotMA(res05, ylim=c(-2,2))

# Volcano Plot
volcanoplot <- ggplot(data = res05_DF, aes(x = log2FoldChange, y = -log10(padj), col=DGE )) +
  geom_vline(xintercept = c(-0, 0), col = "gray", linetype = "dashed")+
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = "dashed") + 
  geom_point(size = 2) +
  scale_color_manual(values = c("blue", "black", "red"),
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  coord_cartesian(ylim = c(0, 65), xlim = c(-8, 8)) + 
  labs(color = 'Genes', 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-adj")) + 
  ggtitle('COVID vs healthy patients')  
#saving volacano plot file
pdf(file = "volcanoplot.pdf", width = 10, height = 10) 
volcanoplot
dev.off()

# PlotCounts
# Top 6 genes with lowest padj
head(res05_DF)
par(mfrow = c(2,3))
plotCounts(dds, gene="GENE", intgroup="Sample_characteristics_ch", main = "GENE") #Replace GENE with specific GENE

## -----------------------------------------------------------------------------
sessionInfo()

