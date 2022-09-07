
# **sev - The SummarizedExperiment Viewer**

## Overview
The sev package is a utility package that relies completely on the bioconductor packages *DEP* and *iSEE* / *iSEEu* for the analysis, adding convenience functions, imputation methods, GO-term enrichment via *clusterProfiler* and additional *iSEE* containers.

## Installation
Open up R and run the following code:

```{r eval = FALSE, echo=T}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("DO.db")

#All other dependencies are installed automatically with sev. 
devtools::install_github("adhutz/sev")
```

## View vignette
An examplary analysis is included as a vignette. 
```{r eval = FALSE, echo=T}
vignette("sev")
```

## Run example
The examplary analysis is also included in the installation folder as an Rmarkdown document together with the proteinGroups.txt file. To open the corresponding folder directly, run the following code.
```{r eval = FALSE, echo=T}
utils::browseURL(system.file("extdata", package = "sev"))
```

