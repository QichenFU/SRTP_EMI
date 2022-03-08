genes <- read.table("gene.txt",header = T)
ep_cell <- read.table("all_ep_stem.txt",header = T)
EMI_cell <- read.csv("EMI_stem_cell_colcell.csv",row.names = genes$Index)
EMI_cell <- EMI_cell[,2:ncol(EMI_cell)]
row.names(ep_cell) <- genes$Index

ep_cell_names <- names(ep_cell)
EMI_cell_names <- names(EMI_cell)

ep_cell_not_normal <- ep_cell[,!grepl("LUNG_N",ep_cell_names)]
EMI_cell_not_normal <- EMI_cell[,!grepl("LUNG_N",EMI_cell_names)]

library(Seurat)
# clustering
ep_object <- CreateSeuratObject(counts = ep_cell_not_normal)
ep_object <- FindVariableFeatures(ep_object, selection.method = "vst", nfeatures = 2000)
ep_object <- ScaleData(ep_object)
ep_object <- RunPCA(ep_object)
DimHeatmap(ep_object, dims = 1:6, cells = 500, balanced = TRUE)

ep_object <- JackStraw(ep_object, num.replicate = 100)
ep_object <- ScoreJackStraw(ep_object, dims = 1:20)
JackStrawPlot(ep_object, dims = 1:15)
ElbowPlot(ep_object)

ep_object <- FindNeighbors(ep_object, dims = 1:10)
ep_object <- FindClusters(ep_object, resolution = 0.5)

ep_object <- RunUMAP(ep_object, dims = 1:10)
DimPlot(ep_object, reduction = "umap",label = TRUE)

umap_embeddings <- data.frame(ep_object@reductions[["umap"]]@cell.embeddings)
temp <- c()
for (i in row.names(umap_embeddings)){
  if (i %in% names(ep_cell_not_normal)){
    temp <- c(temp,1)
  }
  else{
    temp <- c(temp,0)
  }
}
umap_embeddings$label <- temp

library(ggplot2)
g <- ggplot(umap_embeddings,aes(x=UMAP_1,y=UMAP_2,color=label))
g <- g + geom_point()
g

umap_embeddings$VIM <- t(ep_cell_not_normal["VIM",])[,1]
umap_embeddings$CDH1 <- t(ep_cell_not_normal["CDH1",])[,1]
g <- ggplot(umap_embeddings,aes(x=UMAP_1,y=UMAP_2,color=VIM))
g <- g + geom_point()+ scale_color_gradient(low = "cyan",high = "red")
g
g <- ggplot(umap_embeddings,aes(x=UMAP_1,y=UMAP_2,color=CDH1))
g <- g + geom_point() + scale_color_gradient(low = "cyan",high = "red")
g
g <- ggplot(umap_embeddings,aes(x=UMAP_1,y=UMAP_2,color=CDH1))
g <- g + geom_point(aes(size=VIM))
g <- g + scale_color_gradientn(colours = c('cyan','red'))
g
library(Nebulosa)
plot_density(ep_object,c("VIM", "CDH1"))
# https://baijiahao.baidu.com/s?id=1714013229702372282&wfr=spider&for=pc
library(dplyr)

ep_object.markers <- FindAllMarkers(ep_object,logfc.threshold=0.25)
ep_object.markers_for_geo <- ep_object.markers[which(ep_object.markers$cluster==1&ep_object.markers$p_val_adj<0.05&ep_object.markers$avg_log2FC>1.5),"gene"]
# cluster 0: esophagus (S100A2)
# cluster 1: 
# cluster 3: lung (SCGB3A1, MEG3)

ep_object.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top6 <- ep_object.markers %>% group_by(cluster) %>% top_n(n = 6, wt = avg_log2FC)
DoHeatmap(ep_object, features = top6$gene) + NoLegend()
bottom6 <- ep_object.markers %>% group_by(cluster) %>% top_n(n = 6, wt = -avg_log2FC) 
DoHeatmap(ep_object, features = bottom6$gene) + NoLegend()


