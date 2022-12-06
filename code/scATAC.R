library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(patchwork)
set.seed(1234)
#aPSM.counts <- Read10X_h5("/Volumes/macbook_backup/backup/star_protocol/aPSM_scATAC/filtered_peak_bc_matrix.h5")
#aPSM_meta <- read.table("/Volumes/macbook_backup/backup/star_protocol/aPSM_scATAC/singlecell.csv.gz", sep = ",", header = TRUE, row.names = 1)
aPSM.counts <- Read10X_h5("./aPSM_scATAC/filtered_peak_bc_matrix.h5")
aPSM_meta <- read.table("./aPSM_scATAC/singlecell.csv.gz", sep = ",", header = TRUE, row.names = 1)
aPSM_chrom_assay <- CreateChromatinAssay(
  counts = aPSM.counts,
  sep = c(":","-"),
  genome = 'mm10',
  fragments = './aPSM_scATAC/fragments.tsv.gz', min.cells = 3, min.features = 100)
aPSM_atac <- CreateSeuratObject(
  counts = aPSM_chrom_assay, 
  assay = 'aPSM_peaks',
  project = 'aPSM_atac',
  meta.data = aPSM_meta[colnames(aPSM_chrom_assay),])

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"
# add the gene information to the object
Annotation(aPSM_atac) <- annotations

aPSM_atac <- NucleosomeSignal(object = aPSM_atac)
# compute TSS enrichment score per cell
aPSM_atac <- TSSEnrichment(object = aPSM_atac, fast = FALSE)
aPSM_atac$pct_reads_in_peaks <- aPSM_atac$peak_region_fragments / aPSM_atac$passed_filters * 100
aPSM_atac$blacklist_ratio <- aPSM_atac$blacklist_region_fragments / aPSM_atac$peak_region_fragments
#tiff("../figures/figure8.tiff", units="in", width=10, height=10, res=300)
#VlnPlot(
#  object = aPSM_atac,
#  features = c('pct_reads_in_peaks', 'peak_region_fragments',
#               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),pt.size = 0.1, ncol = 3)
#dev.off()

#jpeg("../figures/figure8.jpg", units="in", width=10, height=5, res=300)
VlnPlot(
  object = aPSM_atac,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'nucleosome_signal'),pt.size = 0.1,ncol=4)
#dev.off()



FeatureScatter(aPSM_atac, feature1 = "peak_region_fragments", feature2 = "nCount_aPSM_peaks")

aPSM_atac <- subset(
  x = aPSM_atac,
  subset = peak_region_fragments > 2586 &
    peak_region_fragments < 20000 & pct_reads_in_peaks > 15 & blacklist_ratio < 0.05)
ncol(aPSM_atac)
VlnPlot(
  object = aPSM_atac,
  features = c('nucleosome_signal','peak_region_fragments'),pt.size = 0.1) + NoLegend()

aPSM_atac <- RunTFIDF(aPSM_atac)
aPSM_atac <- FindTopFeatures(aPSM_atac, min.cutoff = 'q0')
aPSM_atac <- RunSVD(
  object = aPSM_atac, assay = 'aPSM_peaks',
  reduction.key = 'LSI_', reduction.name = 'lsi' )

library(ggplot2)
aPSM_atac <- RunUMAP(object = aPSM_atac, reduction = 'lsi', dims = 1:40)
aPSM_atac <- RunTSNE(object = aPSM_atac, reduction = 'lsi', dims = 1:40)
aPSM_atac <- FindNeighbors(object = aPSM_atac, reduction = 'lsi', dims = 1:40)

aPSM_atac <- FindClusters(object = aPSM_atac, verbose = FALSE,resolution=0.25)
tiff("../figures/figure9.tiff", units="in", width=10, height=5, res=300)
DimPlot(object = aPSM_atac, label = F,reduction = 'umap') +labs(title = " aPSM scATAC")
dev.off()

aPSM_gene.activities <- GeneActivity(aPSM_atac)
save(aPSM_gene.activities,file="aPSM_atac_gene.activities.RData")
# add the gene activity matrix to the Seurat object as a new assay and normalize it
aPSM_atac[['RNA']] <- CreateAssayObject(counts = aPSM_gene.activities)
aPSM_atac <- NormalizeData(
  object = aPSM_atac, assay = 'RNA', normalization.method = 'LogNormalize',
  scale.factor = median(aPSM_atac$nCount_RNA) )
save(aPSM_atac,file="aPSM_scATAC.RData")





