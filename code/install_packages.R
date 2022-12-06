# Package names
packages <- c("ggplot2","dplyr","patchwork","Rcurl")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")
BiocManager::install("Rsamtools")
BiocManager::install("EnsDb.Mmusculus.v79")
BiocManager::install("biovizBase")
BiocManager::install("chromVAR")
BiocManager::install("motifmatchr")
BiocManager::install("JASPAR2020")
BiocManager::install("TFBSTools")
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")

install.packages("R.utils")
install.packages("remotes")

## install SeuratWrappers
remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)

