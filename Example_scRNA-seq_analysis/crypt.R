library(dplyr)
library(Seurat)
library(patchwork)
library(scater)
library(SingleCellExperiment)
library(pheatmap)

wtintest.data <- Read10X(data.dir = "D:/Interns/2020.9-/sc RNA/data and code/crypt_normal")
wtintest <- CreateSeuratObject(counts = wtintest.data, project = "wtintest",min.cells = 3,min.features = 200)

injuriedintest.data <- Read10X(data.dir = "D:/Interns/2020.9-/sc RNA/data and code/crypt_irradiated")
injuriedintest <- CreateSeuratObject(counts = injuriedintest.data, project = "injuriedintest",,min.cells = 3,min.features = 200)

# crypt_normal_data <- Read10X(data.dir = "D:/Interns/2020.9-/sc RNA/data and code/crypt_normal")
# crypt_normal <- CreateSeuratObject(counts = crypt_normal_data, project = "crypt_normal")

merged <- merge(wtintest,y = injuriedintest,add.cell.ids = c("wt","injuried"),project = "merged")

sce <- SingleCellExperiment(crypt_normal)
sce

# mito DNA 
wtintest.data[c("mt-Nd1",
                "mt-Nd2",
                "mt-Co1",
                "mt-Co2",
                "mt-Atp8",
                "mt-Atp6",
                "mt-Co3",
                "mt-Nd3",
                "mt-Nd4l",
                "mt-Nd4",
                "mt-Nd5",
                "mt-Nd6",
                "mt-Cytb"
),1]
dense.size <- object.size(as.matrix(wtintest.data))
dense.size
wtintest[["percent.mt"]] <- PercentageFeatureSet(wtintest, pattern = "mt-")
injuriedintest[["percent.mt"]] <- PercentageFeatureSet(injuriedintest, pattern = "mt-")
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "mt-")

VlnPlot(wtintest, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(injuriedintest, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(wtintest, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(wtintest, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

wtintest <- subset(wtintest, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
injuriedintest <- subset(injuriedintest, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
merged <- subset(merged, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)

wtintest <- NormalizeData(wtintest, normalization.method = "LogNormalize", scale.factor = 10000)
injuriedintest <- NormalizeData(injuriedintest, normalization.method = "LogNormalize", scale.factor = 10000)
#wtintest <- NormalizeData(wtintest)

all.genes <- rownames(merged)
merged <- ScaleData(merged, features = all.genes)
merged <- FindVariableFeatures(object = merged)

wtintest <- FindVariableFeatures(wtintest, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(wtintest), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(wtintest)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(wtintest)
wtintest <- ScaleData(wtintest, features = all.genes)

wtintest <- RunPCA(wtintest, features = VariableFeatures(object = wtintest))
merged <- RunPCA(merged, npcs = 50, features = VariableFeatures(object = merged))
# Examine and visualize PCA results a few different ways
print(wtintest[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(wtintest, dims = 1:2, reduction = "pca")
DimPlot(wtintest, reduction = "pca")
DimHeatmap(wtintest, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(wtintest, dims = 1:15, cells = 500, balanced = TRUE)

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
wtintest <- JackStraw(wtintest, num.replicate = 100)
wtintest <- ScoreJackStraw(wtintest, dims = 1:20)

JackStrawPlot(wtintest, dims = 1:15)

ElbowPlot(wtintest)

wtintest <- FindNeighbors(wtintest, dims = 1:10)
wtintest <- FindClusters(wtintest, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(wtintest), 5)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
wtintest <- RunUMAP(wtintest, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(wtintest, reduction = "umap")
ElbowPlot(wtintest)

# Cluster the cells
wtintest <- FindNeighbors(wtintest, dims = 1:10)
merged <- FindNeighbors(merged, dims = 1:10)
wtintest <- FindClusters(wtintest, resolution = 0.5)
merged <- FindClusters(merged, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(wtintest), 5)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
wtintest <- RunUMAP(wtintest, dims = 1:10)
wtintest <- RunTSNE(wtintest, dims = 1:10)
merged <- RunTSNE(merged,dims=1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(wtintest, reduction = "tsne")

saveRDS(wtintest, file = "D:/Interns/2020.9-/sc RNA/data and code/wtintest3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/wtintest_tutorial.rds")

# find all markers of cluster 1
cluster1.markers <- FindMarkers(wtintest, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(wtintest, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
wtintest.markers <- FindAllMarkers(wtintest, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
merged.markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pheatmap(merged.markers)

wtintest.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

cluster1.markers <- FindMarkers(wtintest, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(wtintest, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
VlnPlot(wtintest, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(wtintest, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                               "CD8A"))

top10 <- wtintest.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(wtintest, features = top10) + NoLegend()
DoHeatmap(merged, features = merged.markers$gene) + NoLegend()

new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(wtintest)
wtintest <- RenameIdents(wtintest, new.cluster.ids)
DimPlot(wtintest, reduction = "tsne", label = TRUE, pt.size = 1) 
DimPlot(merged, reduction = "tsne", label = TRUE, pt.size = 1) 
# saveRDS(wtintest, file = "D:/Interns/2020.9-/sc RNA/data and code/wtintest3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/wtintest3k_final.rds")
