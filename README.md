---
title: "README.md"
author: "Arvid Hutzfeldt"
date: '2022-07-13'
output:
  html_document:
    theme: cosmo
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

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

## Data import
Data can be imported directly from MaxQuant derived proteinGroups.txt or similar files to create a summarized experiment. 