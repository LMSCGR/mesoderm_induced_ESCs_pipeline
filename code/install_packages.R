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

install.packages("R.utils")
install.packages("remotes")

## install SeuratWrappers
remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)

