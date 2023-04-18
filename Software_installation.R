install.packages('Seurat')

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("XVector","Rhtslib","GenomeInfoDb", "GenomicRanges", "IRanges", "Rsamtools", "S4Vectors", "BiocGenerics"))

install.packages('Signac')
install.packages('harmony')
install.packages('Cairo')

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats','limma', 'lme4', 'S4Vectors', 'SingleCellExperiment', 'SummarizedExperiment', 'batchelor','HDF5Array', 'terra', 'ggrastr'))

install.packages('grr')
url="https://cran.r-project.org/src/contrib/Archive/Matrix.utils/"
filename="Matrix.utils_0.9.8.tar.gz"

download.file(paste(url, filename, sep = ""), paste(getwd(),filename,sep = "/"), mode = "wb")

install.packages(paste(getwd(),filename,sep = "/"), repos = NULL, type="source")
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/monocle3')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
BiocManager::install("JASPAR2020")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
BiocManager::install("TFBSTools")

# Package names
packages <- c("ggplot2","dplyr","patchwork","RCurl","readr","grid","gridExtra","hdf5r")
# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {install.packages(packages[!installed_packages])}
BiocManager::install("biomaRt")
BiocManager::install("Rsamtools")
BiocManager::install("EnsDb.Mmusculus.v79")
BiocManager::install("biovizBase")
BiocManager::install("chromVAR")
BiocManager::install("motifmatchr")
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
install.packages("R.utils")
install.packages("remotes")
# install SeuratWrappers
remotes::install_github('satijalab/seurat-wrappers')
