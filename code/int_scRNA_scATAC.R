library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(patchwork)
library(ggplot2)
set.seed(1234)

load("aPSM_scRNA.RData")
load("aPSM_scATAC.RData")

DefaultAssay(aPSM_atac) <- 'RNA'
ncol(aPSM_atac)
transfer.anchors <- FindTransferAnchors(
  reference = aPSM, query = aPSM_atac, k.anchor = 20,
  k.filter = 200, reduction = 'cca', dims = 1:30)
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = aPSM$seurat_clusters,
  weight.reduction = aPSM_atac[['lsi']],
  dims = 2:30)
save(transfer.anchors,file="transfer.anchors_aPSM_atac.RData")
aPSM_atac <- AddMetaData(object = aPSM_atac, metadata = predicted.labels)
save(aPSM_atac,file="aPSM_atac_meta.RData")
DimPlot(object = aPSM_atac, label = F,reduction = 'umap',group.by ='predicted.id' ) +labs(title = " aPSM scATAC")
DimPlot(object = aPSM, label = F,reduction = 'umap') +labs(title = " aPSM")
aPSM_gene<-read.csv("aPSM_f.txt",header =F)[,1]
names(aPSM_gene)<-aPSM_gene
del_list<-c('Rora','Tbx18','Vtn','Sox6','Ripply2','Pgm5','Psck5','Myocd','Met','Fhf18','Dmrt2','Cer1','Foxp1',
            'Zic5','Sim1','Pcsk5','Pappa','Meox2','Gadd45b','Fgf18','Cd36','Rock2','Rock1','Jup','Fyn','Calm1',
            'Actb','Apoe')
aPSM_gene<-aPSM_gene[setdiff(aPSM_gene,del_list )]
aPSM_gene<-unname(aPSM_gene)
p1<-DotPlot(aPSM_atac, features = aPSM_gene, cols = c("blue", "red"), dot.scale = 8, assay="RNA",scale.min =0, scale.max = 100,group.by ='predicted.id' ) +coord_flip() + xlab('Gene') +  ylab('Cluster')+ guides(color = guide_colorbar(title = 'Average estimated expression(z-score)'))+labs(title = "aPSM scATAC")
s1<-DotPlot(aPSM, features = aPSM_gene, cols = c("blue", "red"), dot.scale = 8, assay="RNA", scale.min =0,scale.max = 100) +coord_flip() + xlab('Gene') +  ylab('Cluster')+ guides(color = guide_colorbar(title = 'Average expression(z-score)'))+labs(title = "aPSM scRNA")

library(grid)
library(gridExtra)
tiff("../figures/figure13.tiff", units="in", width=15, height=7.5, res=300)
grid.arrange(s1, p1, ncol = 2,widths = c(8,9))
dev.off()

