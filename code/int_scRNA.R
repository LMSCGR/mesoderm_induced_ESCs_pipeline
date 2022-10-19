library(dplyr)
library(Seurat)
library(harmony)
library(data.table)
library(parallel)
# Set number of cores to use
NCORES = 1 #
meta <- fread("naive_instructed_esc.csv")
data_dir <- list("/Volumes/macbook_backup/backup/star_protocol/naive_scRNA/","/Volumes/macbook_backup/backup/star_protocol/instructed_scRNA/")
mat.list <- list()
soupx.used <- list()
for(i in 1:length(data_dir)){ 
  mat.list[[i]] <- Read10X(data.dir = paste0(data_dir[i], 'filtered_feature_bc_matrix'))
  soupx.used[[i]] <- F}
cat(sum(unlist(lapply(mat.list, ncol))),"cells (total) loaded...\n")

sample_num<-min(ncol(mat.list[[1]]),ncol(mat.list[[2]]))
sel.id<-sample(colnames(mat.list[[2]]), size=sample_num, replace=FALSE)
mat.list[[2]]<-mat.list[[2]][,sel.id]

seu.list <- list()
seu.list <- mclapply( 
  mat.list,
  FUN = function(mat){
    return(CreateSeuratObject(
      counts = mat, min.features = 200, min.cells = 3,project = 'naive_instructed_data'))
  }, mc.cores = NCORES)  
for(i in 1:length(seu.list)){
  cat(' ------------------------------------\n',
      '--- Processing dataset number ', i, '-\n',
      '------------------------------------\n')
  # Add meta data
  for(md in colnames(meta)){
    seu.list[[i]][[md]] <- meta[[md]][i]
  }
  
  # add %MT
  seu.list[[i]][["percent.mt"]]  <- PercentageFeatureSet(seu.list[[i]], pattern = "mt-") 
  
  # Filter out low quality cells according to the metrics defined above
  seu.list[[i]] <- subset(seu.list[[i]],
                          subset = nFeature_RNA > 1600 & nFeature_RNA < 8000 & percent.mt < 20)
  # Only mito and floor filtering; trying to find doublets
}
cat((sum(unlist(lapply(mat.list, ncol)))-sum(unlist(lapply(seu.list, ncol)))),"cells (total) removed...\n")
seuPreProcess <- function(seu, assay='RNA', n.pcs=30, res=0.25){
  
  pca.name = paste0('pca_', assay)
  pca.key = paste0(pca.name,'_')
  umap.name = paste0('umap_', assay)
  
  seu = NormalizeData(
    seu
  ) %>% FindVariableFeatures(
    assay = assay,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = F
  ) %>% ScaleData(
    assay = assay
  ) %>% RunPCA(
    assay = assay,
    reduction.name = pca.name,
    reduction.key = pca.key,
    verbose = F,
    npcs = n.pcs
  )
  n.pcs.use =n.pcs
  # FindNeighbors %>% RunUMAP, FindClusters
  seu <- FindNeighbors(
    seu,
    reduction = pca.name,
    dims = 1:n.pcs.use,
    force.recalc = TRUE,
    verbose = FALSE
  ) %>% RunUMAP(
    reduction = pca.name,
    dims = 1:n.pcs.use,
    reduction.name=umap.name
  )
  seu@reductions[[umap.name]]@misc$n.pcs.used <- n.pcs.use  
  seu <- FindClusters(object = seu,resolution = res)
  seu[[paste0('RNA_res.',res)]] <- as.numeric(seu@active.ident)
  return(seu)
}
# preprocess each dataset individually
seu.list <- lapply(seu.list, seuPreProcess)

tmp.list <- list()
for(i in 1:length(seu.list)){
  DefaultAssay(seu.list[[i]]) <- "RNA"
  tmp.list[[i]] <- DietSeurat(seu.list[[i]], assays = "RNA")
}

# merge tmp count matrices
scMuscle.pref.seurat <- merge(
  tmp.list[[1]],
  y = tmp.list[[2]]
)
tiff("../figures/figure10.tiff", units="in", width=10, height=5, res=300)
VlnPlot(
  scMuscle.pref.seurat,
  features = c(
    'nCount_RNA',
    'nFeature_RNA',
    'percent.mt'
  ),
  group.by = 'source',
  pt.size = 0
)
dev.off()

#     Seurat preprocessing on merged data ----
DefaultAssay(scMuscle.pref.seurat) <- 'RNA'
scMuscle.pref.seurat <- 
  NormalizeData(
    scMuscle.pref.seurat, assay = 'RNA'
  ) %>% FindVariableFeatures(
    selection.method = 'vst',
    nfeatures = 2000,
    verbose = TRUE
  ) %>% ScaleData(
    assay = 'RNA',
    verbose = TRUE
  ) %>% RunPCA(
    assay = 'RNA',
    reduction.name = 'pca_RNA',
    reduction.key = 'pca_RNA_',
    verbose = TRUE,
    npcs = 50
  ) 
tiff("../figures/figure11.tiff", units="in", width=10, height=5, res=300)
ElbowPlot(scMuscle.pref.seurat, reduction = 'pca_RNA', ndims = 50)
dev.off()

n.pcs = 30
scMuscle.pref.seurat <-
  RunUMAP(
    scMuscle.pref.seurat, reduction = 'pca_RNA', 
    dims = 1:n.pcs, reduction.name='umap_RNA'
  ) %>% FindNeighbors(
    reduction = 'pca_RNA',
    dims = 1:n.pcs,
    force.recalc = TRUE,
    verbose = F
  ) 

scMuscle.pref.seurat <- FindClusters(object = scMuscle.pref.seurat, resolution = 0.25)
scMuscle.pref.seurat[['RNA_res.0.25']] <- as.numeric(scMuscle.pref.seurat@active.ident)

scMuscle.pref.seurat <- 
  scMuscle.pref.seurat %>% RunHarmony(
    group.by.vars=c('sample'), reduction='pca_RNA',
    assay='RNA',plot_convergence = TRUE,verbose=TRUE) 
n.pcs = npcs(scMuscle.pref.seurat,reduction="harmony")
scMuscle.pref.seurat <- 
  scMuscle.pref.seurat %>% RunUMAP(
    reduction = 'harmony', dims = 1:n.pcs,
    reduction.name='umap_harmony')  
scMuscle.pref.seurat@reductions$umap_harmony@misc$n.pcs.used <- n.pcs
scMuscle.pref.seurat <- 
  scMuscle.pref.seurat %>% FindNeighbors(
    reduction = 'harmony',dims = 1:n.pcs, 
    graph.name = 'harmony_snn',force.recalc = TRUE,
    verbose = FALSE)
scMuscle.pref.seurat <- FindClusters(
  object = scMuscle.pref.seurat,resolution = 1.0,
  graph.name='harmony_snn')
scMuscle.pref.seurat[['harmony_res.1.0']] <- as.numeric(scMuscle.pref.seurat@active.ident) 
scMuscle.pref.seurat <- FindClusters(
  object = scMuscle.pref.seurat,
  resolution = 2.0, graph.name='harmony_snn')
scMuscle.pref.seurat[['harmony_res.2.0']] <- as.numeric(scMuscle.pref.seurat@active.ident)

library(cowplot)
library(ggplot2)

p1<-DimPlot(object = scMuscle.pref.seurat, reduction = "umap_RNA", pt.size = .1, group.by = "sample")+labs(title = "Merged by Seurat directly")
p2<-DimPlot(object = scMuscle.pref.seurat, reduction = "umap_harmony", pt.size = .1, group.by = "sample")+labs(title = "Merged by Seurat with Harmony")
tiff("../figures/figure12.tiff", units="in", width=10, height=5, res=300)
p1+p2
dev.off()
save(scMuscle.pref.seurat,file="naive_instructed_scRNA_ESCs.RData")




