# Install required packages (if not already installed)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("RankProd", quietly = TRUE)) {
  BiocManager::install("RankProd")
}
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  BiocManager::install("tidyverse")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  BiocManager::install("dplyr")
}
if (!requireNamespace("magrittr", quietly = TRUE)) {
  BiocManager::install("magrittr")
}
if (!requireNamespace("stringr", quietly = TRUE)) {
  BiocManager::install("stringr")
}
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

# Load packages
library(RankProd)
library(tidyverse)
library(dplyr)
library(magrittr)
library(stringr)
library(pheatmap)
library(RColorBrewer)

## -----------------------------------------------------------------------------
# Import results from previous DGE analyses
res1 <- read.csv("res05_1.csv")
res2 <- read.csv("res05_2.csv")
res3 <- read.csv("res05_3.csv")
res4 <- read.csv("res05_4.csv")

# (If required) Process res separately due to duplicates
res4 <- aggregate(res4, FUN = mean, by = c(list(res4$GeneName), list(res4$DGE)))
res4 <- res4[,-c(3,10)]
colnames(res4) <- c("GeneName", "DGE", "baseMean", "log2FoldChange", 
                      "lfcSE", "stat", "pvalue", "padj")

# Obtain LFC from results
LFC1 <- res1[,c(1,3)]
LFC2 <- res2[,c(1,3)]
LFC3 <- res3[,c(1,3)]
LFC4 <- res4[,c(1,4)]

# Merge all LGC data
LFC <- LFC1 %>%
  full_join(LFC2, by= "GeneName") %>%
  full_join(LFC3, by= "GeneName") %>%
  full_join(LFC4, by= "GeneName") 
colnames(LFC) <- c("GeneName", "Study1", "Study2", "Study3", "Study4")
row.names(LFC) <- LFC[,1]
LFC <- LFC[,-1]
write.table(LFC, "swab_LFC.txt")

## -----------------------------------------------------------------------------
# Filter NA values
na_max <- apply(as.matrix(LFC), 1, function(x) sum(is.na(x)))
LFC_f <- LFC[na_max <= 3,]

# Functions
# Replace NA with mean 
na2average <- function(x) {
    if (any(is.na(x))) {
        x[is.na(x)] <- mean(x, na.rm = TRUE)
    }
    x
}


## -----------------------------------------------------------------------------
# Rank Calculations
# Calculate ranks for downregulated genes
ranks.down.mx <- apply(LFC_f, 2, rank, ties.method="average", na.last="keep")
ranks.down <- t(apply(ranks.down.mx, 1, na2average))
ranks.down.mx2 <- apply(ranks.down, 2, rank, ties.method = "average", na.last = "keep")

# Calculating ranks for upregulated genes
ranks.up.mx <- apply((LFC_f * -1), 2, rank, ties.method="average", na.last="keep")
ranks.up <- t(apply(ranks.up.mx, 1, na2average))
ranks.up.mx2 <- apply(ranks.up, 2, rank, ties.method = "average", na.last = "keep")

# Plot correlations
# Function for Pearson correlation plotting
panel.cor <- function(x, y, digits = 3, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, use="pairwise.complete.obs")
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0("r=", txt)
  # if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  # text(0.5, 0.5, txt, cex = cex.cor * r)
  text(0.5, 0.5, txt, cex=1.5, col="blue", font = 3)
}

# Function for Spearman correlation plotting
panel.cor.spearman <- function(x, y, digits = 3, prefix = "", cex.cor, ...) {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, use = "pairwise.complete.obs", method = "spearman")
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0("r=", txt)
    text(0.5, 0.5, txt, cex = 1.5, col = "purple", font = 3)
}

# Execute Pearson correlation plot
cor.df<- LFC %>% 
  set_colnames(names(LFC)) 
colnames(cor.df) <- str_replace(colnames(cor.df), "Study1", "GSE152075")
colnames(cor.df) <- str_replace(colnames(cor.df), "Study2", "GSE156063")
colnames(cor.df) <- str_replace(colnames(cor.df), "Study3", "GSE163151")
colnames(cor.df) <- str_replace(colnames(cor.df), "Study4", "GSE166530") 

pairs(cor.df, col = "brown", pch=20, upper.panel = panel.cor,
      xaxt='n', yaxt='n',
      text.panel = function(x,y,lab,cex,font) {text(x,y,lab, cex=1.0, font=2)}) 

# Execute Spearman correlation plot
ranks.up.mx2 %>% 
  set_colnames(names(LFC)) 
colnames(ranks.up.mx2) <- str_replace(colnames(ranks.up.mx2), "Study1", "GSE152075")
colnames(ranks.up.mx2) <- str_replace(colnames(ranks.up.mx2), "Study2", "GSE156063")
colnames(ranks.up.mx2) <- str_replace(colnames(ranks.up.mx2), "Study3", "GSE163151")
colnames(ranks.up.mx2) <- str_replace(colnames(ranks.up.mx2), "Study4", "GSE166530") 

pairs(ranks.up.mx2, col = "brown", pch=18, upper.panel = panel.cor.spearman,
      xaxt='n', yaxt='n',
      text.panel = function(x,y,lab,cex,font) {text(x,y,lab, cex=1.5, font=2)})

ranks.down.mx2 %>% 
  set_colnames(names(LFC)) 
colnames(ranks.down.mx2) <- str_replace(colnames(ranks.down.mx2), "Study1", "GSE152075")
colnames(ranks.down.mx2) <- str_replace(colnames(ranks.down.mx2), "Study2", "GSE156063")
colnames(ranks.down.mx2) <- str_replace(colnames(ranks.down.mx2), "Study3", "GSE163151")
colnames(ranks.down.mx2) <- str_replace(colnames(ranks.down.mx2), "Study4", "GSE166530") 

pairs(ranks.down.mx2, col = "brown", pch=18, upper.panel = panel.cor.spearman,
      xaxt='n', yaxt='n',
      text.panel = function(x,y,lab,cex,font) {text(x,y,lab, cex=1.5, font=2)})

## -----------------------------------------------------------------------------
# RankProd Analysis Functions
# Function to calculate rankprod
get.rank.prod <- function(mx) {
  cl.down <- rep(1, ncol(mx))
  rank.prod <- RankProd::RankProducts(
    data=mx,
    cl=cl.down,
    calculateProduct=TRUE,
    gene.names = row.names(ranks.down.mx2)
  )
  rank.prod
}

get.rank.prod.stats <- function(rank.prod.result) {
  # Get rank prods:
  RP <- as.data.frame(rank.prod.result$RPs) %>% 
    tibble::rownames_to_column(.) %>% 
    dplyr::select(1:2) %>% 
    set_colnames(c("Gene", "Rank.Prod")) %>% 
    arrange(Rank.Prod)
  
  # Get p values:
  PVAL <- as.data.frame(rank.prod.result$pval) %>% 
    tibble::rownames_to_column(.) %>% 
    dplyr::select(1:2) %>% 
    set_colnames(c("Gene", "P.Val"))
  
  # Merge:
  RP.PVAL <- merge(x=RP, y=PVAL, by="Gene", all=T)
  
  # Adjust p values:
  RP.PVAL$P.Val.Adj <- p.adjust(RP.PVAL$P.Val, method = "fdr")
  
  RP.PVAL
}

## -----------------------------------------------------------------------------
# RankProd Analysis Execution
# Calculate rank products for up and down regulated genes
rank.prod.up <- get.rank.prod(ranks.up.mx2) %>% 
  get.rank.prod.stats(.)
rank.prod.down <- get.rank.prod(ranks.down.mx2) %>% 
  get.rank.prod.stats(.)

# Get top 10 genes
top.10.up <- rank.prod.up %>% 
  arrange(Rank.Prod) %>% 
  head(., 10) 
top.10.down <- rank.prod.down %>% 
  arrange(Rank.Prod) %>% 
  head(., 10) 

# Save results
write.table(top.10.up, file = "top10_up_regulated.txt", sep = "\t")
write.table(top.10.down, file = "top10_down_regulated.txt", sep = "\t")

## -----------------------------------------------------------------------------
# Process Significant Genes
# Get number of significant down-regulated genes
n_sig_down <- nrow(subset(rank.prod.down, P.Val.Adj <= 0.05))
# Get number of significant up-regulated genes
n_sig_up <- nrow(subset(rank.prod.up, P.Val.Adj <= 0.05))

# Get top 20 down-regulated genes
down_top20 <- dplyr::arrange(rank.prod.down, P.Val.Adj, Rank.Prod) %>% 
    head(20) %>% 
    .$Gene
mx_down_top20 <- ranks.down.mx2[down_top20,]

# Get top 20 up-regulated genes
up_top20 <- dplyr::arrange(rank.prod.up, P.Val.Adj, Rank.Prod) %>% 
    head(20) %>% 
    .$Gene
mx_up_top20 <- ranks.up.mx2[up_top20,]

# Save Significant Genes
# Process and save significantly up-regulated genes
sig_up_genes <- rank.prod.up %>% 
    dplyr::filter(P.Val.Adj <= 0.05) %>%
    arrange(P.Val.Adj)
write.table(sig_up_genes, 
            file = "significant_upregulated_genes.txt", 
            sep = "\t", 
            row.names = FALSE)
# Process and save significantly down-regulated genes
sig_down_genes <- rank.prod.down %>% 
    dplyr::filter(P.Val.Adj <= 0.05) %>%
    arrange(P.Val.Adj)
write.table(sig_down_genes, 
            file = "significant_downregulated_genes.txt", 
            sep = "\t", 
            row.names = FALSE)

## -----------------------------------------------------------------------------
# Visualisation
# Generate Color Scheme
heat_colors <- colorRampPalette(RColorBrewer::brewer.pal(7, "Blues"))(256) %>% rev

# Generate and Save Heatmaps
pdf("differential_expression_heatmaps.pdf", width = 12, height = 8)

# Up-regulated genes heatmap
pheatmap(mx_up_top20[1:10,],
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 10,
         main = "Top 10 Up-regulated Genes",
         color = heat_colors)

# Down-regulated genes heatmap
pheatmap(mx_down_top20[1:10,],
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 10,
         main = "Top 10 Down-regulated Genes",
         color = heat_colors)

dev.off()

sessionInfo()
