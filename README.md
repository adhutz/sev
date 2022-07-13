
# **sev - The SummarizedExperiment Viewer**

## Overview
The sev package is a utility package that relies completely on the bioconductor packages *DEP* and *iSEE* / *iSEEu* for the analysis, adding convenience functions, imputation methods, GO-term enrichment via *clusterProfiler* and additional *iSEE* containers.

## Installation
Open up R and run the following code:

```{r eval = FALSE, echo=T}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("iSEE")
BiocManager::install("iSEEu")

#All other dependencies are installed automatically with sev. 
devtools::install_github("adhutz/sev")
```


