library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
set.seed(1234)

wt_q_60h_sce<-as.SingleCellExperiment(wt_q_60h)
library(slingshot)
colData(wt_q_60h_sce)
wt_q_60h_sce <- slingshot(wt_q_60h_sce, clusterLabels = 'cell_type',reducedDim = "PCA",
                      allow.breaks = FALSE)
library(gam)
t <- order(wt_q_60h_sce$slingPseudotime_1)

# for time, only look at the 100 most variable genes 
Y <- log1p(assay(wt_q_60h_sce,"logcounts"))

var100 <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:100]
Y <- Y[var100,]

# fit a GAM with a loess term for pseudotime
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  suppressWarnings({
    tmp <- gam(z ~ lo(t), data=d)
  })
  p <- summary(tmp)[3][[1]][2,3]
  p
})

heatdata <- assays(wt_q_60h_sce)$logcounts[b_genelist, order(t, na.last = NA)]
heatclus <- wt_q_60h_sce$cell_type[order(t, na.last = NA)]
library(Polychrome)
my_color <- createPalette(10, c("#010101", "#ff0000"), M=1000)
names(my_color) <- unique(as.character(wt_q_60h_sce$cell_type))

heatmap(as.matrix(heatdata), Colv = NA,
        ColSideColors = my_color[heatclus],cexRow = 1,cexCol = 1)

#####monocle
library(monocle)
load("/Users/kok3/Desktop/development/proj_ezh1/figures/data/qui_60h_p_monocle_final_new.RData")
q_60h<-q_60h_p_samplesheet[q_60h_p_samplesheet$timepoint!="P",]
cell_id<-rownames(q_60h)
b_genelist <-
  row.names(subset(fData(q_60h_p_CellData_pseudo),
                   gene_short_name %in% c("Ezh1","Ezh2","Myog")))

plot_pseudotime_heatmap(q_60h_p_CellData_pseudo[b_genelist,cell_id],
                        num_clusters = 2,
                        norm_method ="vstExprs",
                        cores = 1,
                        show_rownames = T)

####seurat
q_60h_exp<-exprs(q_60h_p_CellData_pseudo)[,cell_id]

wt_q_60h <- CreateSeuratObject(counts = q_60h_exp, project = "wt_q_60h")

wt_q_60h  <- NormalizeData(object = wt_q_60h, normalization.method = "LogNormalize", scale.factor = 1e4)
max(rowMeans(as.matrix(wt_q_60h@assays$RNA@data)))
wt_q_60h <- FindVariableFeatures(wt_q_60h, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
wt_q_60h_top10 <- head(VariableFeatures(wt_q_60h), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(wt_q_60h)
plot2 <- LabelPoints(plot = plot1, points = wt_q_60h_top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

wt_q_60h.all.genes <- rownames(wt_q_60h)
wt_q_60h <- ScaleData(wt_q_60h, features = wt_q_60h.all.genes)

wt_q_60h <- RunPCA(object = wt_q_60h, features = VariableFeatures(object = wt_q_60h), verbose = FALSE)
ElbowPlot(object = wt_q_60h,ndims =50)

wt_q_60h <- FindNeighbors(object = wt_q_60h, dims = 1:30)
wt_q_60h <- FindClusters(object = wt_q_60h, resolution = 0.25)
wt_q_60h <- RunTSNE(object = wt_q_60h, dims = 1:30)
wt_q_60h <- RunUMAP(object = wt_q_60h, dims = 1:30)
DimPlot(object = wt_q_60h, reduction = 'umap',label=T)+labs(title = "wt_q_60h")

cell_type<-c(rep("qui",length(wt_q.id)),rep("36h",length(wt_36h.id)),rep("60h",length(wt_60h.id)))
names(cell_type)<-c(wt_q.id,wt_36h.id,wt_60h.id)
cell_type<-factor(cell_type)

wt_q_60h@meta.data$cell_type<-q_60h$timepoint
DimPlot(object = wt_q_60h, reduction = 'umap',group.by='cell_type',label=T)+labs(title = "wt_q_60h")
FeaturePlot(wt_q_60h,features = c("ENSMUSG00000029687", "ENSMUSG00000009471" ,"ENSMUSG00000006920","ENSMUSG00000026459"))
save(wt_q_60h,file="wt_q_60h_seur.RData")
