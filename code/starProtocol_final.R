
rm(list = ls())

library(Seurat)
library(Signac)
library(patchwork)
library(monocle3) 
library(SeuratWrappers)

library(EnsDb.Mmusculus.v79)
library(GenomeInfoDb)

library(ggplot2)
library( dplyr)

set.seed(1234)

# setup working directory
setwd("path/starProtocol")

# Star Protocol 2: Multiomics analysis 

# Data preprocessing
# Load snRNA and snATAC data
Star.data <- Read10X (data.dir = "path/filtered_feature_bc_matrix")

# Extract RNA and ATAC data
rna_counts <- Star.data$`Gene Expression` 
atac_counts <- Star.data$Peaks

# create a Seurat object containing the RNA data
Star <- CreateSeuratObject(
  counts = rna_counts,
  project = "Star", min.cells=5, min.features = 100,
  assay = "RNA")

# Select peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]

# Add annotation
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "mm10"

# Create ATAC assay
fragpath <- "path/Star_atac_fragments.tsv.gz"

# create ATAC assay and add it to the object
Star[["ATAC"]] <- CreateChromatinAssay(
  counts  = atac_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = fragpath,
  min.cells = 10,
  annotation = annotation
)

# Downsize the dataset 
set.seed(111)
Star <- subset(x = Star, downsample = 6000)

# Quality control plots(Figure 5A) for snRNA and snATAC before removing low quality cells
DefaultAssay(Star) <- "RNA"
Star[["percent.mt"]] <- PercentageFeatureSet(Star, pattern = "^mt-")

VlnPlot(Star, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()

DefaultAssay(Star) <- "ATAC"
Star <- NucleosomeSignal(Star)
Star <- TSSEnrichment(object=Star, fast=FALSE)

VlnPlot(Star, features = c("nCount_ATAC", "nFeature_ATAC", 
                           "TSS.enrichment", "nucleosome_signal"), 
        ncol = 4, log = TRUE, pt.size = 0) + NoLegend()

# Low quality cells were removed (Figure 5B).
Star <- subset(
  x = Star,
  subset = nCount_ATAC < 1e5 &
    nCount_ATAC > 1e2 &
    nCount_RNA < 100000 &
    nCount_RNA > 1200 &
    nucleosome_signal < 2.5 &
    TSS.enrichment > 3 &
    percent.mt < 10 )

saveRDS(Star, file="Star.rds")

VlnPlot(Star, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()
VlnPlot(Star, features = c("nCount_ATAC", "nFeature_ATAC", 
                           "TSS.enrichment", "nucleosome_signal"), ncol = 4,
        log = TRUE, pt.size = 0) + NoLegend()

# Perform normalization and dimensional reduction of snRNA-seq and 
# snATAC-seq assays independently and individually.

# snRNA analysis
DefaultAssay(Star) <- "RNA"
Star <- SCTransform(Star, verbose = FALSE) %>% RunPCA() %>%
  RunUMAP(dims = 1:30, reduction.name = 'umap', reduction.key = 'UMAP_')

# snATAC analysis
# The first dimension was excluded, as it is typically correlated with sequencing depth
DefaultAssay(Star) <- "ATAC"
Star <- RunTFIDF(Star)
Star <- FindTopFeatures(Star, min.cutoff = 'q0')
Star <- RunSVD(Star)
Star <- RunUMAP(Star, reduction = 'lsi', dims = 2:30,
                reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# WNN analysis
Star <- FindMultiModalNeighbors(Star, reduction.list = list("pca", "lsi"),
                                dims.list = list(1:30, 2:30))
Star <- RunUMAP(Star, nn.name = "weighted.nn", reduction.name = "umap.wnn",
                reduction.key = "wnnUMAP_")
Star <- FindClusters(Star, graph.name = "wsnn",  resolution = 0.8, algorithm = 3, verbose = FALSE)

# Visualize clustering based on snRNA-seq, snATAC-seq, or WNN analysis (Fig6A).
p1 <- DimPlot(Star, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 8, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(Star, reduction = "umap.atac",group.by = "seurat_clusters", label = TRUE, label.size = 8, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(Star, reduction = "umap.wnn", group.by = "seurat_clusters", label = TRUE, label.size = 8, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() 


# Find markers for every cluster compared to all remaining cells
# Assign the cell type annotation. 
DefaultAssay(Star) <- "RNA"

Star.rna.markers <- FindAllMarkers(Star, assay = "RNA", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Star.rna.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)


# Cell state annotations were added
Star <- RenameIdents(Star, '10' = 'cell_3','11' = 'cell_1','12' = 'cell_2','13' = 'cell_5')
Star <- RenameIdents(Star, '5' = 'cell_4','6' = 'cell_3','7' = 'cell_2','8' = 'cell_5','9' = 'cell_3')
Star <- RenameIdents(Star, '0' = 'cell_1','1' = 'cell_1','2' = 'cell_1','3' = 'cell_2','4' = 'cell_1')

Star$celltype <- Idents(Star)

# visualize clustering based on gene expression, ATAC-seq, or WNN analysis (Figure 6D).
p1 <- DimPlot(Star, reduction = "umap", group.by = "celltype", label = FALSE, label.size = 8, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(Star, reduction = "umap.atac",group.by = "celltype", label = FALSE, label.size = 8, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(Star, reduction = "umap.wnn", group.by = "celltype", label = FALSE, label.size = 8, repel = TRUE) + ggtitle("WNN")

p1+p2+p3

# snATAC-seq analysis
# load library
library(chromVAR)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)

DefaultAssay(Star) <- "ATAC"

Star.atac.markers <- FindAllMarkers(Star, assay = "ATAC", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Star.atac.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)


# Adding motif information to the Seurat object
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))

# Add motif information
Star <- AddMotifs(
  object = Star,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pwm_set, 
  assay="ATAC"
)

# Computing motif activities
Star <- RunChromVAR(
  object = Star,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

saveRDS(Star, file="Star_final.rds")

# Star Protocol 2: Data visualization and interpretation 
Star <- readRDS("Star_final.rds")

DefaultAssay(Star) <- 'RNA'

set.seed(22)
cds <- SeuratWrappers::as.cell_data_set(Star, assay = "RNA", reduction = "umap", group.by = "seurat_clusters")
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(Star[["RNA"]])

cds <- preprocess_cds(cds, method = "PCA")
cds <- reduce_dimension(cds, preprocess_method = "PCA",umap.n_neighbors= 14L,
                        reduction_method = "UMAP")
cds <- cluster_cells(cds, reduction_method = "UMAP")
plot_cells(cds, color_cells_by = "cluster", group_label_size = 4, cell_size = 1)
cds <- learn_graph(cds, use_partition = FALSE, close_loop = FALSE)

cell_ids <- colnames(cds)[Star$seurat_clusters ==  "0"]
closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
closest_vertex <- closest_vertex[cell_ids, ]
closest_vertex <- as.numeric(names(which.max(table(closest_vertex))))
mst <- principal_graph(cds)$UMAP
root_pr_nodes <- igraph::V(mst)$name[closest_vertex]

rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name
cds <- order_cells(cds, root_pr_nodes = root_pr_nodes)

# Visualization: Trajectory plot (Figure 6B)
plot_cells(cds, color_cells_by = "pseudotime",
           label_cell_groups =T, label_leaves = F,
           label_branch_points = F,show_trajectory_graph = T,
           graph_label_size = 3,label_groups_by_cluster = T)

# Visualization: Cell states derived from trajectory inference (Figure 6C)
plot_cells(cds, color_cells_by = "cluster", cell_size = 1,
           label_cell_groups = TRUE,group_label_size = 4,
           show_trajectory_graph = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE)

# Visualization: Paired-plots expression of Myod1 and Myog (Figure 7A)
Star.seur <- as.Seurat(cds, assay = NULL, clusters = "UMAP")
Star.seur<-AddMetaData(Star.seur,metadata= cds@principal_graph_aux$UMAP$pseudotime,
                       col.name = "monocle3_pseudotime")

FeaturePlot(Star.seur,features = c("Myod1","Myog"),
            reduction ="UMAP",combine = T,
            blend = TRUE,
            blend.threshold = 0.0, min.cutoff = 0,
            max.cutoff = 6)

# Footprinting visualization: Motifs of interest were footprinted with default 
# parameters between the cell types. (Figure 7B). 

Star_135 <- subset(x = Star, idents = c("cell_1",  "cell_3", "cell_5"), invert = FALSE)
DefaultAssay(Star_135) <- "ATAC"

# Gather the footprinting information for sets of motifs
Star_135 <- Footprint(
  object = Star_135,
  motif.name = c("MYOG", "MYOD1"),
  genome = BSgenome.Mmusculus.UCSC.mm10)

PlotFootprint(Star_135, features = c("MYOD1")) + patchwork::plot_layout(ncol = 1)
PlotFootprint(Star_135, features = c("MYOG")) + patchwork::plot_layout(ncol = 1)

sessionInfo()
