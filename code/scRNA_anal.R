library(dplyr)
library(Seurat)
library(monocle3)
library(SeuratWrappers)
# for plotting
library(ggplot2)
library(patchwork)
set.seed(1234)
aPSM.matrix <- Read10X(data.dir ="/Volumes/macbook_backup/backup/star_protocol/aPSM_scRNA/filtered_feature_bc_matrix/")
aPSM <- CreateSeuratObject(counts = aPSM.matrix, min.cells = 3, min.features = 200, project = " aPSM")
aPSM[["percent.mt"]] <- PercentageFeatureSet(aPSM, pattern = "^mt-")
jpeg("../figures/figure2a.jpg", units="in", width=10, height=5, res=300)
VlnPlot(aPSM, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
plot1 <- FeatureScatter(aPSM, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(aPSM, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
tiff("../figures/figure4.tiff", units="in", width=10, height=5, res=300)
CombinePlots(plots = list(plot1, plot2))
dev.off()
aPSM<- subset(aPSM, subset = nFeature_RNA > 0 & nFeature_RNA < 8000 & percent.mt < 20)

aPSM <- NormalizeData(object = aPSM, normalization.method = "LogNormalize", scale.factor = 1e4)
aPSM <- FindVariableFeatures(aPSM, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
aPSM_top10 <- head(VariableFeatures(aPSM), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(aPSM)
plot2 <- LabelPoints(plot = plot1, points = aPSM_top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

aPSM.all.genes <- rownames(aPSM)
aPSM <- ScaleData(aPSM, features = aPSM.all.genes)

convertHumanGeneList <- function(x){
  require("biomaRt")
  human 	<- useMart("ensembl", dataset = "hsapiens_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
  mouse 	<- useMart("ensembl", dataset = "mmusculus_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
  tmp 	<- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=TRUE)
  mousex 	<- unique(tmp[,2])
  return(mousex)}
s.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes)
g2m.genes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes)
cell_cycle<-t(read.csv(file="cell_cycle.txt",header=F))[,1]
filtered_genes<-c(s.genes,cell_cycle)


aPSM <- RunPCA(object = aPSM, features = VariableFeatures(object = aPSM), verbose = FALSE)
aPSM <- CellCycleScoring(aPSM, s.features = filtered_genes, g2m.features = g2m.genes, set.ident = TRUE)
aPSM <- ScaleData(aPSM, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(aPSM))
aPSM <- JackStraw(aPSM, num.replicate = 100)
aPSM <- ScoreJackStraw(aPSM, dims = 1:20)
JackStrawPlot(aPSM, dims = 1:20)
tiff("../figures/figure5.tiff", units="in", width=6, height=4, res=300)
ElbowPlot(object = aPSM,ndims =50)
dev.off()

aPSM <- FindNeighbors(object = aPSM, dims = 1:30)
aPSM <- FindClusters(object = aPSM, resolution = 0.25)
aPSM <- RunTSNE(object = aPSM, dims = 1:30)
aPSM <- RunUMAP(object = aPSM, dims = 1:30)
tiff("../figures/figure6.tiff", units="in", width=10, height=5, res=300)
DimPlot(object=aPSM,reduction='umap',label=T)+labs(title = " aPSM")
dev.off()
save(aPSM,file="aPSM_scRNA.RData")

aPSM.markers <- FindAllMarkers(aPSM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
aPSM.markers_table<-aPSM.markers %>%group_by(cluster) %>% slice_max(n = 20, order_by = avg_log2FC)

DefaultAssay(aPSM) <- "RNA"
aPSM_cds<-as.cell_data_set(aPSM)
aPSM_cds<-cluster_cells(aPSM_cds,reduction="UMAP",k = 30,resolution = 0.00012)
aPSM_cds <- learn_graph(aPSM_cds, close_loop = F,use_partition = T,learn_graph_control =list(minimal_branch_len=5))
plot_cells(aPSM_cds, label_groups_by_cluster = T, label_leaves = F, label_branch_points = T,graph_label_size = 3)
aPSM.min.umap <- which.min(unlist(FetchData(aPSM, "UMAP_2")))
aPSM.min.umap <- colnames(aPSM)[aPSM.min.umap]
aPSM_cds <- order_cells(aPSM_cds, root_cells = aPSM.min.umap)
tiff("../figures/figure7.tiff", units="in", width=5, height=3, res=300)
plot_cells(aPSM_cds, color_cells_by = "pseudotime", label_cell_groups =T, label_leaves = F, label_branch_points = F,show_trajectory_graph = T,graph_label_size = 3,label_groups_by_cluster = T)
dev.off()


