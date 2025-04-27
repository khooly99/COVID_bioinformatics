# COVID_bioinformatics
We used a standardised bioinformatics pipeline to run a meta-analysis to identify differentially expressed genes across multiple RNA-Seq datasets of COVID-19 against healthy controls which are grouped into their respective samples. This repository consists R source codes for steps of DGE analysis with DESeq2, meta-analysis with RankProd and WGCNA.

## Dataset
We used publicly available gene expression dataset from Gene Expression Omnibus (GEO) which are as follows:
1. [GSE152075](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152075)
2. [GSE156063](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156063)
3. [GSE163151](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163151)
4. [GSE166530](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166530)
6. [GSE152418](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152418)
7. [GSE179627](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179627)
8. [GSE196822](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196822)
9. [GSE190496](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190496)
10. [GSE202182](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202182)
11. [GSE208076](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208076)
    
## Usage:
1. Install R and R studio:
Install R and R studio from [posit.co](https://posit.co/download/rstudio-desktop/).
2. Install required main packages:
```
install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("RankProd")
install.packages("WGCNA")
```
3. Run following analysis
  - DGE analysis with DESeq2 of individual RNA-Seq datasets from GEO, NCBI.
  - Meta-analysis with RankProd of grouped samples (of swab, blood, and tissue samples)
  - Co-expression analysis with WGCNA of grouped samples (of swab, blood, and tissue samples)
