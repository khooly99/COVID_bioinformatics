# Install required packages (if not already installed)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("WGCNA", quietly = TRUE)) {
  BiocManager::install("WGCNA")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  BiocManager::install("ggplot2")
}
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  BiocManager::install("tidyverse")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  BiocManager::install("dplyr")
}

# Load packages
library(WGCNA)
require(parallel)
require(doParallel)
enableWGCNAThreads()
library(tidyverse)
library(dplyr)
library(impute)
library(limma)
library(ggplot2)

## -----------------------------------------------------------------------------
# Import normalised counts and metadata
# Read normalised counts for datasets
n1 <- read.table("n1.txt", header=TRUE, row.names = NULL)
n2 <- read.table("n2.txt", header=TRUE, row.names = NULL)
n3 <- read.table("n3.txt", header=TRUE, row.names = NULL)
n4 <- read.table("n4.txt", header=TRUE, row.names = NULL)

# Read metadata
n1_meta <- as.data.frame(t(read.csv("n1_meta.csv", header=TRUE, row.names = 1)))
n2_meta <- as.data.frame(t(read.csv("n2_meta.csv", header=TRUE, row.names = 1)))
n3_meta <- as.data.frame(t(read.csv("n3_meta.csv", header=TRUE, row.names = 1)))
n4_meta <- as.data.frame(t(read.csv("n4_meta.csv", header=TRUE, row.names = 1)))

## -----------------------------------------------------------------------------
# Pre-processing of datasets (If required, clean as needed to standardise)
n1_meta$Sample_characteristics_ch[which(n1_meta$Sample_characteristics_ch == "sars-cov-2 positivity: pos")] <- "positive"
n1_meta$Sample_characteristics_ch[which(n1_meta$Sample_characteristics_ch == "sars-cov-2 positivity: neg")] <- "negative"

n2_meta <- n2_meta[match(colnames(n2[,-1]), rownames(n2_meta)),]
n2_meta$Sample_characteristics_ch[which(n2_meta$Sample_characteristics_ch == "sars-cov-2 rpm: POS")] <- "positive"
n2_meta$Sample_characteristics_ch[which(n2_meta$Sample_characteristics_ch == "sars-cov-2 rpm: NEG")] <- "negative"

n3_meta <- n3_meta[match(colnames(n3[,-1]), rownames(n3_meta)),]
n3_meta$Sample_characteristics_ch[which(n3_meta$Sample_characteristics_ch == "pathogen: SARS-CoV-2")] <- "positive"
n3_meta$Sample_characteristics_ch[which(n3_meta$Sample_characteristics_ch == "pathogen: No pathogen")] <- "negative"

n4_meta$Sample_characteristics_ch[which(n4_meta$Sample_characteristics_ch == "subject status: COVID-19 Positive")] <- "positive"
n4_meta$Sample_characteristics_ch[which(n4_meta$Sample_characteristics_ch == "subject status: COVID-19 Negative")] <- "negative"

# Merge normalised counts as grouped samples (In this case, swab)
swabs <- n1 %>%
  full_join(n2, by= "SYMBOL") %>%
  full_join(n3, by= "SYMBOL") %>%
  full_join(n4, by= "SYMBOL")

swabs %>%
  group_by(SYMBOL) %>%
  summarise_all(mean) %>%
  data.frame() -> swabs_df
row.names(swabs_df) <- swabs_df$SYMBOL
swabs_df$SYMBOL <- NULL
swabs_df <- as.data.frame(t(swabs_df))

## -----------------------------------------------------------------------------
# QC
# Missing data removal/imputation
missingValueThreshold <- 0.5
missingProportion <- colMeans(is.na(swabs_df))
missingProportion
swabs.f_df <- swabs_df[,missingProportion <= missingValueThreshold]
swabs.if <- impute.knn(as.matrix(swabs.f_df))$data
swabs.if_df <- as.data.frame(swabs.if)

# Visualisation with PCA before correction
pca <- prcomp(swabs.if_df)
pca.dat <- as.data.frame(pca$x)
pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits=2)
pca <- ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point()+
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1:', pca.var.percent[1], '%'),
       y = paste0('PC2:', pca.var.percent[2], '%'))

# Batch effect correction
# Grouping of batch
sampleID <- as.data.frame(c(row.names(n1_meta), row.names(n2_meta), row.names(n3_meta),
                            row.names(n4_meta)))
colnames(sampleID) <- "Sample"
condition <- as.data.frame(c(n1_meta$Sample_characteristics_ch, n2_meta$Sample_characteristics_ch,
                             n3_meta$Sample_characteristics_ch, n4_meta$Sample_characteristics_ch))
colnames(condition) <- "Condition"
batch <- c(rep(1, ncol(n1[,-1])), rep(2, ncol(n2[,-1])), rep(3, ncol(n3[,-1])), rep(4, ncol(n4[,-1])))
swab_meta <- cbind(sampleID, condition, batch)

# Executing correction
swab_meta$Condition <- as.factor(swab_meta$Condition)
swab_meta$batch <- as.factor(swab_meta$batch)
design <- model.matrix(~0 + swab_meta$Condition)
swabs.cif <- removeBatchEffect(t(swabs.if_df), batch = swab_meta$batch, design = design)
swabs.cif_df <- as.data.frame(t(swabs.cif))

# Check gene outliers with gsg
gsg <- goodSamplesGenes(swabs.cif_df, verbose=3)
gsg$allOK

# Visualisation with PCA after correction
pca.f <- prcomp(swabs.cif_df)
pca.f.dat <- as.data.frame(pca.f$x)
pca.f.var <- pca.f$sdev^2
pca.f.var.percent <- round(pca.f.var/sum(pca.f.var)*100, digits=2)
pca.f <- ggplot(pca.f.dat, aes(PC1, PC2)) +
  geom_point()+
  geom_text(label = rownames(pca.f.dat)) +
  labs(x = paste0('PC1:', pca.f.var.percent[1], '%'),
       y = paste0('PC2:', pca.f.var.percent[2], '%'))

## -----------------------------------------------------------------------------
# WGCNA Network Construction
# Clustering
sampleTree <- hclust(dist(swabs.cif_df), method="average")
par(cex = 0.5)
par(mar = c(2,4,2,2))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h = 300, col = "red")

# Determine clusters below the line
clust <-  cutreeStatic(sampleTree, cutHeight = 300, minSize = 10)
table(clust)
keepSamples <-(clust==1)
swabs.cif_df1 <-  swabs.cif_df[keepSamples,]
nGenes <-  ncol(swabs.cif_df1)
nSamples <-  nrow(swabs.cif_df1)

table(swab_meta$Sample==row.names(swabs.cif_df1))
swab_meta1 <- swab_meta[match(row.names(swabs.cif_df1), swab_meta$Sample),]

# Select soft power threshold
powers <- c(c(1:10), seq(from=10, to=30, by=1)) #choosing a set of soft-thresholding powers
sft <- pickSoftThreshold(swabs.cif_df1, powerVector=powers, verbose=5, networkType="signed") #call network topology analysis function
sft.dat <- sft$fitIndices
par(mfrow= c(1,2))
cex1=0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab= "Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type= "n", main= paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")
abline(h=0.8, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= "Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")

# Build network
softPower <-  15 # Chosen based on the graphs before
temp_cor <- cor
cor <- WGCNA::cor
netwk <- blockwiseModules(swabs.cif_df1,
                           checkMissingData = TRUE,
                           power = softPower,
                           networkType = "signed",
                           minModuleSize = 15,
                           maxBlockSize = 15000,
                           deepSplit = 2,
                           pamRespectsDendro = FALSE,
                           reassignThreshold = 0,
                           mergeCutHeight = 0.15,
                           numericLabels = F,
                           verbose = 3)

## -----------------------------------------------------------------------------
# Check modules from network
table(netwk$colors)
mergedColors <- netwk$colors
table(mergedColors)
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

# Associating modules and phenotypes
# List of modules
module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = mergedColors)
colnames(module_df) <- c("gene", "module")
module_df$module <- paste0("ME", module_df$module)

# Module-Trait Relationship Analysis
# Calculate module eigengenes
MEs <- moduleEigengenes(swabs.cif_df1, mergedColors)$eigengenes
MEs0 <- orderMEs(MEs)
plotEigengeneNetworks(MEs0, "Eigengene dendogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)

# Correlate with traits
row.names(swab_meta1) <- swab_meta1$Sample
traits <-  swab_meta1 %>%
  mutate(disease_state = ifelse (grepl("positive", Condition), 1, 0)) %>%
  dplyr::select(4)
head(traits)
moduleTraitCor <-  cor(MEs0, traits, use= "p")
moduleTraitPvalue <-  corPvalueStudent(moduleTraitCor, nSamples)

# Create correlation heatmap
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(2, 10, 2, 2))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = "COVID-19",
               yLabels = colnames(MEs0),
               ySymbols = colnames(MEs0),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               xLabelsAngle = 0,
               main = paste("Module-trait relationships of swab samples"))

# Choosing modules MEgreen
module_key <- module_df
module_key <- module_key %>% filter(module == "MEgreen")
rownames(module_key)

## -----------------------------------------------------------------------------
# Results analysis
# Student t-test
module_sig <- "MEgreen"
tMEs <- as.data.frame(MEs0[,module_sig])
tMEs$condition <- traits$disease_state
colnames(tMEs) <- c("MEgreen", "condition")
tMEs$condition <- as.factor(tMEs$condition)
ggplot(tMEs, aes(x=condition, y=MEgreen, color=condition)) +
  geom_boxplot()

# Calculate module membership
ModuleMembership <- cor(MEs0, swabs.cif_df1, use = "p")
ModuleMembership.pval <- corPvalueStudent(ModuleMembership, nSamples)

# Calculate gene significance
GeneSignificance <- cor(swabs.cif_df1, traits$disease_state, use = "p")
GS.pval <- corPvalueStudent(GeneSignificance, nSamples)
colnames(GS.pval) <- "GS.pval"
GS.pval %>%
  as.data.frame() %>%
  arrange(GS.pval) %>%
  head(25)

# Selecting genes with MM (>0.8) and GS(>0.2) as threshold
par(mfrow= c(1,1))
verboseScatterplot(abs(ModuleMembership[module_sig, module_key$gene]),
                   abs(GeneSignificance[module_key$gene, 1]),
                   xlab = paste("Module Membership in green module"),
                   ylab = "Gene significance",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "green")
abline(h=0.2, v=0.8, col="red")

## -----------------------------------------------------------------------------
# Obtaining significant genes
# Identify genes with high Module Membership (MM > 0.8)
MMsig <- data.frame(
  abs(ModuleMembership[module_sig, module_key$gene])>0.8, 
  ModuleMembership[module_sig, module_key$gene])
colnames(MMsig) <- c("MMsig", "MM")
MMsig <- MMsig %>%
  dplyr::filter(MMsig == "TRUE")

# Identify genes with high Gene Significance (GS > 0.2)
GSsig <- data.frame(
  abs(GeneSignificance[module_key$gene, 1])>0.2, 
  GeneSignificance[module_key$gene, 1])
colnames(GSsig) <- c("GSsig", "GS")
GSsig <- GSsig %>%
  dplyr::filter(GSsig == "TRUE")

# Combine MM and GS criteria to identify and extract hub genes
sigHubGenes <- data.frame(abs(ModuleMembership[module_sig, module_key$gene])>0.8, abs(GeneSignificance[module_key$gene, 1])>0.2,
                          ModuleMembership[module_sig, module_key$gene], GeneSignificance[module_key$gene, 1])
colnames(sigHubGenes) <- c("MMsig", "GSsig", "MM", "GS")
sigHubGenes <- sigHubGenes %>%
  as.data.frame() %>%
  dplyr::filter(MMsig == "TRUE") %>%
  dplyr::filter(GSsig == "TRUE")
HubGenes <- as.data.frame(row.names(sigHubGenes))
colnames(HubGenes) <- "HubGene"
write.csv(sigHubGenes, "HubGenes_swab.all.csv")

# Match hub genes with their gene significance p-values
ind <- match(HubGenes$HubGene, row.names(GS.pval))
GS.ind <- data.frame(row.names(GS.pval)[ind],GS.pval[ind,])

# Save gene significance values for the module
GS_mod <- as.data.frame(GeneSignificance[module_key$gene, 1])
write.csv(GS_mod, "GS_swab.all.csv")

sessionInfo()
