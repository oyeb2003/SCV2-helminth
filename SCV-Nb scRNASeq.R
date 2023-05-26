#Integrating Nipocov2 and controlcov2 libraries

library(Seurat)
library(ggplot2)
library(patchwork)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggsci)
library(hdf5r)
library(openxlsx)
library(readxl)
library(rio)
library(data.table)
library(magrittr)
library(utils)
library(magrittr)
library(SingleR)
library("SingleCellExperiment")

set.seed(42)
nipocov2 <- readRDS(file = "nipocov2.rds")
controlcov2 <- readRDS(file = "controlcov2.rds")

#Try combining the datasets now. Only make one UMAP and then split by treatment for analysis.
cov2.anchors <- FindIntegrationAnchors(object.list = list(nipocov2, controlcov2), dims = 1:30)
cov2 <- IntegrateData(anchorset = cov2.anchors, dims = 1:30)
DefaultAssay(cov2) <- "integrated"
cov2 <- NormalizeData(cov2)
cov2 <- ScaleData(cov2)
cov2 <- NormalizeData(cov2, assay = "RNA")
cov2 <- ScaleData(cov2, assay = "RNA")

saveRDS(cov2, file = "cov2.rds")
cov2 <- readRDS(file = "cov2.rds")

DefaultAssay(cov2) <- "RNA"
cov2 <- FindVariableFeatures(cov2, selection.method = "vst", nfeatures = 5000)
top10RNA <- head(VariableFeatures(cov2), 10)
plot1 <- VariableFeaturePlot(cov2)
plot2 <-LabelPoints(plot = plot1, points = top10RNA, repel = TRUE)
plot2

DefaultAssay(cov2) <- "integrated"
####Clustering

#PCA
cov2 <- RunPCA(cov2, npcs =30, verbose = FALSE)
print(cov2[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(cov2, dims = 1:2, reduction = "pca")
DimPlot(cov2, reduction = "pca")
DimHeatmap(cov2, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(cov2, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(cov2, ndims = 40, reduction = "pca")


#cluster the cells
cov2 <- FindNeighbors(cov2, dims = 1:13)
cov2 <- FindClusters(cov2, resolution = 0.2)

#try .2 - .8 cluster resolutions to see what makes a better umap

#run non-linear dimensional reduction (UMAP/tSNE)
cov2 <- RunUMAP(cov2, dims = 1:13)
DimPlot(cov2, label = TRUE, reduction = "umap")

cov2 <- FindClusters(cov2, resolution = 0.5)
cov2 <- RunUMAP(cov2, dims = 1:13)
DimPlot(cov2, label = TRUE, reduction = "umap", pt.size = 1)

cov2 <- FindClusters(cov2, resolution = 0.7)
cov2 <- RunUMAP(cov2, dims = 1:13)
Idents(cov2) <- "orig.ident"
DimPlot(cov2, label = TRUE, reduction = "umap", pt.size = 0.7)
Idents(cov2) <- "seurat_clusters"

cov2 <- FindClusters(cov2, resolution = 0.7)
cov2 <- RunUMAP(cov2, dims = 1:15)
DimPlot(cov2, label = TRUE, reduction = "umap")

saveRDS(cov2, file = "cov2update.rds")
cov2 <- readRDS(file = "cov2.rds")


##Find Markers
DefaultAssay(cov2) <- "RNA"
Idents(cov2) <- "seurat_clusters"

#RNA
cov2.markers <- FindAllMarkers(cov2, only.pos = T)
cov2.markers <- cov2.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(cov2, cov2.markers$gene, size = 5)+
  scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#6300ae','#F8766D','#FFCC99')), midpoint = 0, guide = "colourbar", aesthetics = "fill") 
cov2markers <- cov2.markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
cov2markers
export(RP.rna.markers, file = "/Users/howardng/documents/Rewilding Project/Block 1 + Block 2/CORRECTED_RNAmarkers_top25.xlsx")


#SingleR to name clusters



library("SingleR")
library("Seurat")


#Loading and Using SingleR builtin Refs 
immgen.se <- celldex::ImmGenData()  
table(immgen.se$label.main)

#Setting Single cell experiment data as table and object
write.table(as.matrix(GetAssay(object = cov2, slot = "counts")), 
            'counts2.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)

GetAssay(object = cov2, slot = "counts")

cov2matrix <- write.table(cov2@assays[["RNA"]]@counts, file='Counts.tsv', quote=FALSE, sep='\t', col.names = TRUE)
cov2matrix <- read.table(file = "Counts.tsv")
head(cov2matrix)
cov2matrixsub <- cov2matrix[1:10,1:10]
pred.cov2matrixsub <- SingleR(test = cov2matrixsub, ref = immgen.se, assay.type.test=1,
                                  labels = immgen.se$label.main)

###Based on video
####results <- SingleR(test =as.SingleCellExperiment(cov2), ref = immgen.se, labels = immgen.se$label.main)

pred.cov2matrixsub
table(pred.cov2matrixsub$labels)


pred.cov2matrix <- SingleR(test = cov2matrix, ref = immgen.se, assay.type.test=1,
                               labels = immgen.se$label.main)
table(pred.cov2matrix$labels)


predfine.cov2matrix <- SingleR(test = cov2matrix, ref = immgen.se, assay.type.test=1,
                                   labels = immgen.se$label.fine)
table(predfine.cov2matrix$labels)
export(table(predfine.cov2matrix$labels), file = "cov2nippoSingleRfine.xlsx")
export(table(pred.cov2matrix$labels), file = "cov2nippoSingleRmain.xlsx")
export(pred.cov2matrix, file = "SingleRmainMatrix.xlsx")
export(predfine.cov2matrix, file = "SingleRfineMatrix.xlsx")
rownames(pred.cov2matrix)
export(rownames(pred.cov2matrix), file = "Rows_mainmatrix.xlsx")
export(rownames(predfine.cov2matrix), file = "Rows_finematrix.xlsx")

cov2 <- AddMetaData(object = cov2, metadata= pred.cov2matrix$labels, col.name = "singleR.main")
cov2 <- AddMetaData(object = cov2, metadata= predfine.cov2matrix$labels, col.name = "singleR.fine")


Idents(cov2) <- "singleR.main"
DimPlot(cov2, label = TRUE, reduction = "umap", pt.size = 0.7, repel = TRUE)
Idents(cov2) <- "singleR.fine"
Idents(cov2) <- "seurat_clusters"
DimPlot(cov2, label = TRUE, reduction = "umap", pt.size = 1, repel = TRUE)



FeaturePlot(cov2, features = "Il33", pt.size = 1, order = TRUE, split.by = "Condition")
FeaturePlot(cov2, features = "S100a4", pt.size = 0.7, order = TRUE, split.by = "Condition")
FeaturePlot(cov2, features = "Ifitm3", pt.size = 0.7, order = TRUE)
FeaturePlot(cov2, features = "Apoe", pt.size = 0.7, order = TRUE)
FeaturePlot(cov2, features = "C1qa", pt.size = 0.7, order = TRUE)
FeaturePlot(cov2, features = "C1qb", pt.size = 0.7, order = TRUE)

FeaturePlot(cov2, features = "Chil3", pt.size = 0.7, order = TRUE)
FeaturePlot(cov2, features = "Plet1", pt.size = 0.7, order = TRUE)
FeaturePlot(cov2, features = "Mrc1", pt.size = 0.7, order = TRUE)
FeaturePlot(cov2, features = "Siglecf", pt.size = 0.7, order = TRUE, split.by = "orig.ident")
FeaturePlot(cov2, features = "Itgax", pt.size = 0.7, order = TRUE, split.by = )
FeaturePlot(cov2, features = "MN985325.1", pt.size = 0.7, order = TRUE)
FeaturePlot(cov2, features = "Sars", pt.size = 0.7, order = TRUE, split.by = "orig.ident")

FeaturePlot(cov2, features = "Gmcsf", pt.size = 0.7, order = TRUE)
gmcsf
mcsf
#Testing singleR fine labels

library(magrittr)
library(ggplot2)
Idents(cov2) <- 'singleR.fine'
umap_fine = cov2@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% cbind(fine = cov2@meta.data$singleR.fine)

ggplot(umap_fine, aes(x=UMAP_1, y=UMAP_2, color=fine)) + geom_point() + 
  scale_color_manual(values=c("DC (DC.103-11B+24+)" = "blue", 
                              "DC (DC.8-4-11B+)" = "blue",
                             "DC (DC.8-)" = "blue",
                              "DC (DC.PDC.8+)" = "blue",
                              "DC (DC.8-4-11B-)" = "blue",
                              "Epithelial cells (Ep.8wk.CEC.Sca1+)" = "green2",
                              "Endothelial cells (BEC)" = "red",
                              "Epithelial cells (MECHI.GFP-.ADULT)" = 
                                "green2",
                             "Epithelial cells (Ep.5wk.MEC.Sca1+)" = "green2",
                             "Epithelial cells (MECHI.GFP+.ADULT)" = "green2"))

ggplot(umap_fine, aes(x=UMAP_1, y=UMAP_2, color=fine)) + geom_point() + 
  scale_color_manual(values=c("T cells (T.CD4.1H)" = "blue", 
                              "T cells (T.CD4.24H)" = "blue",
                              "T cells (T.CD4.48H)" = "blue",
                              "T cells (T.CD4.5H)" = "blue",
                              "T cells (T.CD4.CTR)" = "blue",
                              "T cells (T.CD4+TESTDB)" = "blue",
                              "T cells (T.CD4+TESTNA)" = "blue",
                              "T cells (T.CD4CONTROL)" = "blue",
                              "T cells (T.CD4TESTCJ)" = "red",
                              "T cells (T.CD4TESTJS)" = "blue",
                              "Tgd (Tgd.mat.VG2+)" = "green2"))

ggplot(umap_fine, aes(x=UMAP_1, y=UMAP_2, color=fine)) + geom_point() + 
  scale_color_manual(values=c("Fibroblasts (FI.MTS15+)" = "blue", 
                              "Fibroblasts (FI)" = "blue",
                              "Fibroblasts (FRC)" = "blue",
                              "ILC (ILC2)" = "red",
                              "ILC (LPL.NCR+CNK)" = "red",
                              "Macrophages (MF.103-11B+24-)" = "green2",
                              "Macrophages (MF.11C-11B+)" = "green2",
                              "Macrophages (MF.11CLOSER.SALM3)" = "green2",
                              "Macrophages (MF.480INT.NAIVE)" = "green2",
                              "Macrophages (MF.MEDL)" = "green2",
                              "Macrophages (MF.SBCAPS)" = "green2",
"Macrophages (MF)" = "green2",
"Macrophages (MFAR-)" = "green2",
"Macrophages (MFIO5.II-480INT)" = "green2",
"Macrophages (MFIO5.II+480INT)" = "green2"))

ggplot(umap_fine, aes(x=UMAP_1, y=UMAP_2, color=fine)) + geom_point() + 
  scale_color_manual(values=c(
                              "Stromal cells (ST.31-38-44-)" = "blue",
                              "Neutrophils (GN.URAC)" = "red",
                              "Neutrophils (GN.Thio)" = "red",
                              "Neutrophils (GN.ARTH)" = "red",
                              "Monocytes (MO.6+2+)" = "green2",
                              "Monocytes (MO.6C-II+)" = "green2",
                              "Monocytes (MO.6C+II-)" = "green2",
                              "Monocytes (MO.6C+II+)" = "green2"))

ggplot(umap_fine, aes(x=UMAP_1, y=UMAP_2, color=fine)) + geom_point() + 
  scale_color_manual(values=c(
    "NK cells (NK.DAP10-)" = "blue",
    "T cells (T.8EFF.OT1.48HR.LISOVA)" = "red",
    "T cells (T.8EFF.OT1.D10LIS)" = "red",
    "T cells (T.8EFF.OT1.D5.VSVOVA)" = "red",
    "T cells (T.8EFF.OT1LISO)" = "red",
    "T cells (T.8EFF.TBET-.OT1LISOVA)" = "red",
    "T cells (T.8MEM.OT1.D45.LISOVA)" = "green2",
    "T cells (T.8MEMKLRG1-CD127+.D8.LISOVA)" = "green2",
"T cells (T.CD4TESTCJ)" = "yellow",
"T cells (T.Tregs)" = "orange"))

#Rename cell clusters
Idents(cov2) <- "seurat_clusters"
new.cluster.ids <- c("Inflammatory Mac", "Mono.Mac", "T", "Stromal.Fibroblasts", "Monocytes", "Neutrophils", "Endothelial.1", "Endothelial.2", 
                     "Epithelial", "Pro.T.1", "Pro.T.2", "B", "Alv.Mac", "DC", "Fibroblasts")

Idents(cov2) <- "orig.ident"

# Rename identity classes
cov2 <- RenameIdents(cov2, "nipocov2" = "nippocov2", "controlcov2" = "controlcov2")

names(new.cluster.ids) <- levels(cov2)
cov2 <- RenameIdents(cov2, new.cluster.ids)
DimPlot(cov2, label = TRUE, reduction = "umap", pt.size = 1, repel = TRUE)
cov2[["CellType"]] <- Idents(object = cov2)

Idents(cov2) <- "orig.ident"
new.cluster.ids <- c("nippocov2", "controlcov2")


names(new.cluster.ids) <- levels(cov2)
cov2 <- RenameIdents(cov2, new.cluster.ids)
DimPlot(cov2, label = TRUE, reduction = "umap", pt.size = 1, repel = TRUE)
cov2[["Condition"]] <- Idents(object = cov2)


#viral RNA
#orf1ab, S, orf3a, E, M, orf6, orf7a, orf8, N, orf10

viralgenes <- c("orf1ab", "S", "orf3a", "E", "M", "orf6", "orf7a", "orf8", "N", "orf10")

FeaturePlot(cov2, feature = "Itgax", pt.size = 0.7, order = TRUE)
FeaturePlot(cov2, feature = "S", pt.size = 0.7, order = TRUE)
FeaturePlot(cov2, feature = "Cxcr6", pt.size = 0.7, order = TRUE, split.by = "orig.ident")
FeaturePlot(cov2, feature = "Itga4", pt.size = 2, order = TRUE, split.by = "orig.ident")
FeaturePlot(cov2, feature = "Cxcl16", pt.size = 0.7, order = TRUE, split.by = "orig.ident")
FeaturePlot(cov2, feature = "Ccl9", pt.size = 0.7, order = TRUE, split.by = "orig.ident")
FeaturePlot(cov2, feature = "Ccl5", pt.size = 0.7, order = TRUE, split.by = "orig.ident")
FeaturePlot(cov2, feature = "Ccl3", pt.size = 0.7, order = TRUE, split.by = "orig.ident")
FeaturePlot(cov2, feature = "Ccl7", pt.size = 0.7, order = TRUE, split.by = "orig.ident")
FeaturePlot(cov2, feature = "Cxcl9", pt.size = 0.7, order = TRUE, split.by = "orig.ident")
FeaturePlot(cov2, feature = "Cxcl10", pt.size = 0.7, order = TRUE, split.by = "orig.ident")
FeaturePlot(cov2, feature = "Cx3cl1", pt.size = 0.7, order = TRUE, split.by = "orig.ident")
FeaturePlot(cov2, feature = "Ch25h", pt.size = 0.7, order = TRUE, split.by = "orig.ident")

FeaturePlot(cov2, feature = "E", pt.size = 0.7, order = TRUE)
FeaturePlot(cov2, feature = "M", pt.size = 0.7, order = TRUE)
FeaturePlot(cov2, feature = "orf6", pt.size = 0.7, order = TRUE)
FeaturePlot(cov2, feature = "orf7a", pt.size = 0.7, order = TRUE)
FeaturePlot(cov2, feature = "orf8", pt.size = 0.7, order = TRUE)
FeaturePlot(cov2, feature = "N", pt.size = 0.7, order = TRUE)
FeaturePlot(cov2, feature = "orf10", pt.size = 0.7, order = TRUE)

#AddModuleScore can be used to see if expression of given gene set is 
#enriched vs set of randomly selected (but based on expression bins) control genes. 
#This might help to clean up the plot as it sounds like the enrichment of the whole gene set 
#would likely be cell type specific whereas one particular gene 
#might also be expressed in other cell types.

DefaultAssay(cov2) <- "RNA"
cov2_marker_gene_list <- list(c("orf1ab", "S", "orf3a", "E", "M", "orf6", "orf7a", "orf8", "N", "orf10"))
cov2 <- AddModuleScore(cov2, features = cov2_marker_gene_list, name = "cov2genes_score")
FeaturePlot(object = cov2, features = "cov2genes_score1", pt.size = 1, order = TRUE, split.by = "orig.ident")



#Differential Gene Expression

DefaultAssay(cov2) <- "Intergrated"
Idents(cov2) <- "CellType"
Macrophages <- subset(cov2, idents = "Macrophages")
Idents(Macrophages) <- "Condition"
avg.Macrophages <- log1p(AverageExpression(Macrophages, verbose = FALSE)$RNA)
head(avg.Macrophages, n = 20)
avg.Macrophages <-data.frame(avg.Macrophages)
is.data.frame(avg.Macrophages)

Idents(cov2) <- "CellType"
Mono.Mac <- subset(cov2, idents = "Mono.Mac")
Idents(Mono.Mac) <- "Condition"
avg.Mono.Mac <- log1p(AverageExpression(Mono.Mac, verbose = FALSE)$RNA)
head(avg.Mono.Mac, n = 20)
avg.Mono.Mac <-data.frame(avg.Mono.Mac)
is.data.frame(avg.Mono.Mac)

Idents(cov2) <- "CellType"
T <- subset(cov2, idents = "T")
Idents(T) <- "Condition"
avg.T <- log1p(AverageExpression(T, verbose = FALSE)$RNA)
head(avg.T, n = 20)
avg.T <-data.frame(avg.T)
is.data.frame(avg.T)

Idents(cov2) <- "CellType"
Stromal.Fibroblasts <- subset(cov2, idents = "Stromal.Fibroblasts")
Idents(Stromal.Fibroblasts) <- "Condition"
avg.Stromal.Fibroblasts <- log1p(AverageExpression(Stromal.Fibroblasts, verbose = FALSE)$RNA)
head(avg.Stromal.Fibroblasts, n = 20)
avg.Stromal.Fibroblasts <-data.frame(avg.Stromal.Fibroblasts)
is.data.frame(avg.Stromal.Fibroblasts)

#Plot2
genes.to.label = c("Ly6a", "Ly6c1", "Fabp5", "Alox5ap", "Ms4a4c", "Cldn10", "Gbp2", "Gbp7", "Pdcd1lg2", "Jun")
p5 <- ggplot(avg.Macrophages, aes(nippocov2, controlcov2)) + geom_point(color = "steelblue2") + ggtitle("Macrophages")
p5 <- LabelPoints(plot = p5, points = genes.to.label, repel = TRUE, max.overlaps = Inf)
p6 <- ggplot(avg.Mono.Mac, aes(nippocov2, controlcov2)) + geom_point(color = "steelblue2") + ggtitle("Mono.Mac")
p6 <- LabelPoints(plot = p6, points = genes.to.label, repel = TRUE, max.overlaps = Inf)
p7 <- ggplot(avg.T, aes(nippocov2, controlcov2)) + geom_point(color = "steelblue2") + ggtitle("T")
p7 <- LabelPoints(plot = p7, points = genes.to.label, repel = TRUE, max.overlaps = Inf)
p8 <- ggplot(avg.Stromal.Fibroblasts, aes(nippocov2, controlcov2)) + geom_point(color = "steelblue2") + ggtitle("Stromal.Fibroblasts")
plot_grid(p5, p6, p7, p8)

Idents(cov2) <- "CellType"
Monocytes <- subset(cov2, idents = "Monocytes")
Idents(Monocytes) <- "Condition"
avg.Monocytes <- log1p(AverageExpression(Monocytes, verbose = FALSE)$RNA)
head(avg.Monocytes, n = 20)
avg.Monocytes <-data.frame(avg.Monocytes)
is.data.frame(avg.Monocytes)

Idents(cov2) <- "CellType"
Neutrophils <- subset(cov2, idents = "Neutrophils")
Idents(Neutrophils) <- "Condition"
avg.Neutrophils <- log1p(AverageExpression(Neutrophils, verbose = FALSE)$RNA)
head(avg.Neutrophils, n = 20)
avg.Neutrophils <-data.frame(avg.Neutrophils)
is.data.frame(avg.Neutrophils)

Idents(cov2) <- "CellType"
Endothelial.1 <- subset(cov2, idents = "Endothelial.1")
Idents(Endothelial.1) <- "Condition"
avg.Endothelial.1 <- log1p(AverageExpression(Endothelial.1, verbose = FALSE)$RNA)
head(avg.Endothelial.1, n = 20)
avg.Endothelial.1 <-data.frame(avg.Endothelial.1)
is.data.frame(avg.Endothelial.1)

Idents(cov2) <- "CellType"
Endothelial.2 <- subset(cov2, idents = "Endothelial.2")
Idents(Endothelial.2) <- "Condition"
avg.Endothelial.2 <- log1p(AverageExpression(Endothelial.2, verbose = FALSE)$RNA)
head(avg.Endothelial.2, n = 20)
avg.Endothelial.2 <-data.frame(avg.Endothelial.2)
is.data.frame(avg.Endothelial.2)

#Plot2
genes.to.label = c("Ly6a", "Ly6c1", "Fabp5", "Alox5ap", "Ms4a4c", "Cldn10", "Gbp2", "Gbp7", "Pdcd1lg2", "Jun")
p5 <- ggplot(avg.Monocytes, aes(nippocov2, controlcov2)) + geom_point(color = "steelblue2") + ggtitle("Monocytes")
p5 <- LabelPoints(plot = p5, points = genes.to.label, repel = TRUE, max.overlaps = Inf)
p6 <- ggplot(avg.Neutrophils, aes(nippocov2, controlcov2)) + geom_point(color = "steelblue2") + ggtitle("Neutrophils")
p6 <- LabelPoints(plot = p6, points = genes.to.label, repel = TRUE, max.overlaps = Inf)
p7 <- ggplot(avg.Endothelial.1, aes(nippocov2, controlcov2)) + geom_point(color = "steelblue2") + ggtitle("Endothelial.1")
p7 <- LabelPoints(plot = p7, points = genes.to.label, repel = TRUE, max.overlaps = Inf)
p8 <- ggplot(avg.Endothelial.2, aes(nippocov2, controlcov2)) + geom_point(color = "steelblue2") + ggtitle("Endothelial.2")
plot_grid(p5, p6, p7, p8)

Idents(cov2) <- "CellType"
Epithelial <- subset(cov2, idents = "Epithelial")
Idents(Epithelial) <- "Condition"
avg.Epithelial <- log1p(AverageExpression(Epithelial, verbose = FALSE)$RNA)
head(avg.Epithelial, n = 20)
avg.Epithelial <-data.frame(avg.Epithelial)
is.data.frame(avg.Epithelial)

Idents(cov2) <- "CellType"
Pro.T.1 <- subset(cov2, idents = "Pro.T.1")
Idents(Pro.T.1) <- "Condition"
avg.Pro.T.1 <- log1p(AverageExpression(Pro.T.1, verbose = FALSE)$RNA)
head(avg.Pro.T.1, n = 20)
avg.Pro.T.1 <-data.frame(avg.Pro.T.1)
is.data.frame(avg.Pro.T.1)

Idents(cov2) <- "CellType"
Pro.T.2 <- subset(cov2, idents = "Pro.T.2")
Idents(Pro.T.2) <- "Condition"
avg.Pro.T.2 <- log1p(AverageExpression(Pro.T.2, verbose = FALSE)$RNA)
head(avg.Pro.T.2, n = 20)
avg.Pro.T.2 <-data.frame(avg.Pro.T.2)
is.data.frame(avg.Pro.T.2)

Idents(cov2) <- "CellType"
B <- subset(cov2, idents = "B")
Idents(B) <- "Condition"
avg.B <- log1p(AverageExpression(B, verbose = FALSE)$RNA)
head(avg.B, n = 20)
avg.B <-data.frame(avg.B)
is.data.frame(avg.B)

#Plot2
genes.to.label = c("Ly6a", "Ly6c1", "Fabp5", "Alox5ap", "Ms4a4c", "Cldn10", "Gbp2", "Gbp7", "Pdcd1lg2", "Jun")
p5 <- ggplot(avg.Epithelial, aes(nippocov2, controlcov2)) + geom_point(color = "steelblue2") + ggtitle("Epithelial")
p5 <- LabelPoints(plot = p5, points = genes.to.label, repel = TRUE, max.overlaps = Inf)
p6 <- ggplot(avg.Pro.T.1, aes(nippocov2, controlcov2)) + geom_point(color = "steelblue2") + ggtitle("Pro.T.1")
p6 <- LabelPoints(plot = p6, points = genes.to.label, repel = TRUE, max.overlaps = Inf)
p7 <- ggplot(avg.Pro.T.2, aes(nippocov2, controlcov2)) + geom_point(color = "steelblue2") + ggtitle("Pro.T.2")
p7 <- LabelPoints(plot = p7, points = genes.to.label, repel = TRUE, max.overlaps = Inf)
p8 <- ggplot(avg.B, aes(nippocov2, controlcov2)) + geom_point(color = "steelblue2") + ggtitle("B")
plot_grid(p5, p6, p7, p8)


Idents(cov2) <- "CellType"
Transitional.Mac <- subset(cov2, idents = "Transitional.Mac")
Idents(Transitional.Mac) <- "Condition"
avg.Transitional.Mac <- log1p(AverageExpression(Transitional.Mac, verbose = FALSE)$RNA)
head(avg.Transitional.Mac, n = 20)
avg.Transitional.Mac <-data.frame(avg.Transitional.Mac)
is.data.frame(avg.Transitional.Mac)

Idents(cov2) <- "CellType"
DC <- subset(cov2, idents = "DC")
Idents(DC) <- "Condition"
avg.DC <- log1p(AverageExpression(DC, verbose = FALSE)$RNA)
head(avg.DC, n = 20)
avg.DC <-data.frame(avg.DC)
is.data.frame(avg.DC)

Idents(cov2) <- "CellType"
Fibroblasts <- subset(cov2, idents = "Fibroblasts")
Idents(Fibroblasts) <- "Condition"
avg.Fibroblasts <- log1p(AverageExpression(Fibroblasts, verbose = FALSE)$RNA)
head(avg.Fibroblasts, n = 20)
avg.Fibroblasts <-data.frame(avg.Fibroblasts)
is.data.frame(avg.Fibroblasts)

Idents(cov2) <- "CellType"
B <- subset(cov2, idents = "B")
Idents(B) <- "Condition"
avg.B <- log1p(AverageExpression(B, verbose = FALSE)$RNA)
head(avg.B, n = 20)
avg.B <-data.frame(avg.B)
is.data.frame(avg.B)

#Plot2
genes.to.label = c("Ly6a", "Ly6c1", "Fabp5", "Alox5ap", "Ms4a4c", "Cldn10", "Gbp2", "Gbp7", "Pdcd1lg2", "Jun")
p5 <- ggplot(avg.Transitional.Mac, aes(nippocov2, controlcov2)) + geom_point(color = "steelblue2") + ggtitle("Transitional.Mac")
p5 <- LabelPoints(plot = p5, points = genes.to.label, repel = TRUE, max.overlaps = Inf)
p6 <- ggplot(avg.DC, aes(nippocov2, controlcov2)) + geom_point(color = "steelblue2") + ggtitle("DC")
p6 <- LabelPoints(plot = p6, points = genes.to.label, repel = TRUE, max.overlaps = Inf)
p7 <- ggplot(avg.Fibroblasts, aes(nippocov2, controlcov2)) + geom_point(color = "steelblue2") + ggtitle("Fibroblasts")
p7 <- LabelPoints(plot = p7, points = genes.to.label, repel = TRUE, max.overlaps = Inf)
p8 <- ggplot(avg.B, aes(nippocov2, controlcov2)) + geom_point(color = "steelblue2") + ggtitle("B")
plot_grid(p5, p6, p7)

Idents(cov2) <- "CellType"
C0C1dge <- FindMarkers(cov2, ident.1 = "0", ident.2 = "1", verbose = FALSE)
head(environ.B.Fo, n = 50)
table <- head(environ.B.Fo, n = 50)
export(table, file = "/Users/howardng/Documents/Environ.B.Fo Top 50 Genes.xlsx")
rownames(head(environ.B.Fo, n = 50))


#DGE of chow vs hp in Cluster 1, 2, 6, 9, 10, 14, 15, 16
Idents(cov2) <- "seurat_clusters"
zero <- subset(cov2, idents = "0")
Idents(zero) <- "Condition"
avg.zero <- log1p(AverageExpression(zero, verbose = FALSE)$RNA)
head(avg.zero, n = 20)
avg.zero <-data.frame(avg.zero)
is.data.frame(avg.zero)
zero.de.markers <- FindMarkers(zero, ident.1 = "nippocov2", ident.2 = "controlcov2")
head(zero.de.markers, n = 300)
export(zero.de.markers, file = "Cluster0markersnippovscontrol.xlsx")

Idents(cov2) <- "seurat_clusters"
one <- subset(cov2, idents = "1")
Idents(one) <- "Condition"
one.de.markers <- FindMarkers(one, ident.1 = "nippocov2", ident.2 = "controlcov2")
head(one.de.markers, n = 300)

Idents(cov2) <- "seurat_clusters"
two <- subset(cov2, idents = "2")
Idents(two) <- "Condition"
two.de.markers <- FindMarkers(two, ident.1 = "nippocov2", ident.2 = "controlcov2")
head(two.de.markers, n = 300)

Idents(cov2) <- "seurat_clusters"
three <- subset(cov2, idents = "3")
Idents(three) <- "Condition"
three.de.markers <- FindMarkers(three, ident.1 = "nippocov2", ident.2 = "controlcov2")
head(three.de.markers, n = 300)

Idents(cov2) <- "seurat_clusters"
four <- subset(cov2, idents = "4")
Idents(four) <- "Condition"
four.de.markers <- FindMarkers(four, ident.1 = "nippocov2", ident.2 = "controlcov2")
head(four.de.markers, n = 300)

Idents(cov2) <- "seurat_clusters"
five <- subset(cov2, idents = "5")
Idents(five) <- "Condition"
five.de.markers <- FindMarkers(five, ident.1 = "nippocov2", ident.2 = "controlcov2")
head(five.de.markers, n = 300)

Idents(cov2) <- "seurat_clusters"
six <- subset(cov2, idents = "6")
Idents(six) <- "Condition"
six.de.markers <- FindMarkers(six, ident.1 = "nippocov2", ident.2 = "controlcov2")
head(six.de.markers, n = 300)

Idents(cov2) <- "seurat_clusters"
seven <- subset(cov2, idents = "7")
Idents(seven) <- "Condition"
seven.de.markers <- FindMarkers(seven, ident.1 = "nippocov2", ident.2 = "controlcov2")
head(seven.de.markers, n = 300)

Idents(cov2) <- "seurat_clusters"
eight <- subset(cov2, idents = "8")
Idents(eight) <- "Condition"
eight.de.markers <- FindMarkers(eight, ident.1 = "nippocov2", ident.2 = "controlcov2")
head(eight.de.markers, n = 300)

Idents(cov2) <- "seurat_clusters"
nine <- subset(cov2, idents = "9")
Idents(nine) <- "Condition"
nine.de.markers <- FindMarkers(nine, ident.1 = "nippocov2", ident.2 = "controlcov2")
head(nine.de.markers, n = 300)

Idents(cov2) <- "seurat_clusters"
ten <- subset(cov2, idents = "10")
Idents(ten) <- "Condition"
ten.de.markers <- FindMarkers(ten, ident.1 = "nippocov2", ident.2 = "controlcov2")
head(ten.de.markers, n = 300)

Idents(cov2) <- "seurat_clusters"
eleven <- subset(cov2, idents = "11")
Idents(eleven) <- "Condition"
eleven.de.markers <- FindMarkers(eleven, ident.1 = "nippocov2", ident.2 = "controlcov2")
head(eleven.de.markers, n = 300)

Idents(cov2) <- "seurat_clusters"
twelve <- subset(cov2, idents = "12")
Idents(twelve) <- "Condition"
twelve.de.markers <- FindMarkers(twelve, ident.1 = "nippocov2", ident.2 = "controlcov2")
head(twelve.de.markers, n = 300)

Idents(cov2) <- "seurat_clusters"
thirteen <- subset(cov2, idents = "13")
Idents(thirteen) <- "Condition"
thirteen.de.markers <- FindMarkers(thirteen, ident.1 = "nippocov2", ident.2 = "controlcov2")
head(thirteen.de.markers, n = 300)

Idents(cov2) <- "seurat_clusters"
fourteen <- subset(cov2, idents = "14")
Idents(fourteen) <- "Condition"
fourteen.de.markers <- FindMarkers(fourteen, ident.1 = "nippocov2", ident.2 = "controlcov2")
head(fourteen.de.markers, n = 300)

Idents(cov2) <- "CellType"

library(dittoSeq)

dittoBarPlot(
  object = cov2,
  var = "CellType",
  group.by = "Condition",
  do.hover = T,
  theme =theme_classic( base_size = 15
  ))

# What proportion of cells are in each cluster?
prop.table(table(Idents(cov2)))

# How does cluster membership vary by replicate?
table(Idents(cov2), cov2$Condition)

prop.table <- prop.table(table(Idents(cov2), cov2$Condition), margin = 2)
prop.table

prop.table <-data.frame((prop.table))
export(prop.table, file = "proprotion.COV2.xlsx")

#viral RNA
#orf1ab, S, orf3a, E, M, orf6, orf7a, orf8, N, orf10
DefaultAssay(cov2) <- "RNA"
chemokine_marker_gene_list <- list(c("Cxcl16", "Ccl17", "Ccl24", "Ccl9"))
cov2 <- AddModuleScore(cov2, features = chemokine_marker_gene_list, name = "chemokinegenes_score")
FeaturePlot(object = cov2, features = "chemokinegenes_score1", pt.size = 1, order = TRUE, split.by = "orig.ident")


Idents(cov2) <- "Condition"
VlnPlot(cov2, features = 
          "Cxcl16")

data <- data.frame(FetchData(object = cov2, slot = 'data',  vars = c('Cxcl16', 'Condition', 'CellType'), group.by = "Condition"))
is.data.frame(data)
export(data, file = "expressionCxcl16groupbycondition.xlsx", rowNames = T, colNames = T)

VlnPlot(cov2, features = 
          "Ccl17")
data <- data.frame(FetchData(object = cov2, slot = 'data',  vars = c('Ccl17', 'CellType'), group.by = "CellType"))
is.data.frame(data)
export(data, file = "expressionCcl17.xlsx", rowNames = T, colNames = T)

VlnPlot(cov2, features = 
          "Cxcl16", split.by = "Condition")
Idents(cov2) <- "Condition"
VlnPlot(cov2, features = 
          "chemokinegenes_score1", split.by = "Condition")
Idents(cov2) <- "CellType"
VlnPlot(cov2, features = 
          "Cxcl16", split.by = "Condition")
VlnPlot(cov2, features = 
          "orf3a", split.by = "Condition")
VlnPlot(cov2, features = 
          "E", split.by = "Condition")
VlnPlot(cov2, features = 
          "M", split.by = "Condition")
VlnPlot(cov2, features = 
          "orf6", split.by = "Condition")
VlnPlot(cov2, features = 
          "orf7a", split.by = "Condition")
VlnPlot(cov2, features = 
          "orf8", split.by = "Condition")
VlnPlot(cov2, features = 
          "N", split.by = "Condition")
VlnPlot(cov2, features = 
          "orf10", split.by = "Condition")
Idents(cov2) <- "Condition"
VlnPlot(cov2, features = 
          "orf1ab", split.by = "Condition")
