#General scRNAseq workflow example

#Load the libraries you will need
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
library(SingleR)
library("SingleCellExperiment")
library(dittoSeq)
library(ComplexHeatmap)
library(circlize)


#Integrating Nipocov2 and controlcov2 libraries

nipo <- readRDS(file = "nipo.rds")
control <- readRDS(file = "control.rds")

#Try combining the datasets now. Only make one UMAP and then split by treatment for analysis.
anchors <- FindIntegrationAnchors(object.list = list(nipo, control), dims = 1:30)
data <- IntegrateData(anchorset = anchors, dims = 1:30)
DefaultAssay(data) <- "integrated"
data <- NormalizeData(data)
data <- ScaleData(data)
data <- NormalizeData(data, assay = "RNA")
data <- ScaleData(data, assay = "RNA")

dim(cntl)
table(cntl$orig.ident)

##identify variable features for each dataset set
DefaultAssay(cntl) <- "RNA"
cntl <- FindVariableFeatures(cntl, selection.method = "vst", nfeatures = 5000)
top10RNA <- head(VariableFeatures(cntl), 10)
plot1 <- VariableFeaturePlot(cntl)
plot2 <-LabelPoints(plot = plot1, points = top10RNA, repel = TRUE)
plot2

DefaultAssay(cntl) <- "integrated"
####Clustering

#PCA
cntl <- RunPCA(cntl)
print(cntl[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(cntl, dims = 1:2, reduction = "pca")
DimPlot(cntl, reduction = "pca")
####split.by = "ident"
DimHeatmap(cntl, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(cntl, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(cntl, ndims = 40, reduction = "pca")


#cluster the cells
cntl <- FindNeighbors(cntl, dims = 1:13)
cntl <- FindClusters(cntl, resolution = 0.2)
cntl
#try .2 - .8 cluster resolutions to see what makes a better umap

#run non-linear dimensional reduction (UMAP/tSNE)
cntl <- RunUMAP(cntl, dims = 1:13)
DimPlot(cntl, label = TRUE, reduction = "umap")

cntl <- FindClusters(cntl, resolution = 0.5)
cntl <- RunUMAP(cntl, dims = 1:13)
DimPlot(cntl, label = TRUE, reduction = "umap", pt.size = 1)

cntl <- FindClusters(cntl, resolution = 0.7)
cntl <- RunUMAP(cntl, dims = 1:13)
DimPlot(cntl, label = TRUE, reduction = "umap", pt.size = 0.7)
ggsave("ggplot2save.svg", width = 10, height = 10, dpi=700)

Idents(cntl)
Idents(cntl) <- "CellType"
Idents(cntl) <- "orig.ident"
Idents(cntl) <- "seurat_clusters"

DimPlot(cntl, label = TRUE, reduction = "umap", pt.size = 0.7)
Idents(cntl) <- "seurat_clusters"

cntl <- FindClusters(cntl, resolution = 0.7)
cntl <- RunUMAP(cntl, dims = 1:13)
DimPlot(cntl, label = TRUE, reduction = "umap")


##Find Markers
DefaultAssay(cntl) <- "RNA"
Idents(cntl) <- "seurat_clusters"

#RNA
cntl.markers <- FindAllMarkers(cntl, only.pos = T)
cntl.markers <- cntl.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
DoHeatmap(cntl, cntl.markers$gene, size = 5)+
  scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#6300ae','#F8766D','#FFCC99')), midpoint = 0, guide = "colourbar", aesthetics = "fill") 
cntl.rna.markers <- cntl.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
export(cntl.rna.markers, file = "CORRECTED_RNAmarkers_top25.xlsx")




library("SingleCellExperiment")

Idents(cntl) <- "orig.ident"
Idents(cntl) <- "seurat_clusters"

immgen.se <- ImmGenData()
immgen.se
table(immgen.se$label.main)

write.table(as.matrix(GetAssay(object = cntl, slot = "counts")), 
            'counts2.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)

GetAssay(object = cntl, slot = "counts")

cntlmatrix <- write.table(cntl@assays[["RNA"]]@counts, file='Counts.tsv', quote=FALSE, sep='\t', col.names = TRUE)
cntlmatrix <- read.table(file = "Counts.tsv")
head(cntlmatrix)
cntlmatrixsub <- cntlmatrix[1:10,1:10]
pred.cntlmatrixsub <- SingleR(test = cntlmatrixsub, ref = immgen.se, assay.type.test=1,
                                  labels = immgen.se$label.main)


pred.cntlmatrixsub
table(pred.cntlmatrixsub$labels)


pred.cntlmatrix <- SingleR(test = cntlmatrix, ref = immgen.se, assay.type.test=1,
                               labels = immgen.se$label.main)
table(pred.cntlmatrix$labels)


predfine.cntlmatrix <- SingleR(test = cntlmatrix, ref = immgen.se, assay.type.test=1,
                                   labels = immgen.se$label.fine)
table(predfine.cntlmatrix$labels)
export(table(predfine.cntlmatrix$labels), file = "cntlnippoSingleRfine.xlsx", row.Names = T, col.Names = T, quote = F)
export(table(pred.cntlmatrix$labels), file = "cntlnippoSingleRmain.xlsx", row.Names = T, col.Names = T)
export(pred.cntlmatrix, file = "SingleRmainMatrix.xlsx", row.Names = T, col.Names = T)
export(predfine.cntlmatrix, file = "SingleRfineMatrix.xlsx", row.Names = T, col.Names = T)
rownames(pred.cntlmatrix)
export(rownames(pred.cntlmatrix), file = "Rows_mainmatrix.xlsx", row.Names = T, col.Names = T)
export(rownames(predfine.cntlmatrix), file = "Rows_finematrix.xlsx",row.Names = T, col.Names = T)


cntl <- AddMetaData(object = cntl, metadata= pred.cntlmatrix$labels, col.name = "singleR.main")
cntl <- AddMetaData(object = cntl, metadata= predfine.cntlmatrix$labels, col.name = "singleR.fine")


Idents(cntl) <- "CellType"
Idents(cntl) <- "orig.ident"
Idents(cntl) <- "RNA"
Idents(cntl) <- "seurat_clusters"
Idents(cntl) <- "singleR.main"



DimPlot(cntl, label = TRUE, reduction = "umap", pt.size = 0.7, repel = TRUE, split.by = "orig.ident")
Idents(cntl) <- "singleR.main"
Idents(cntl) <- "seurat_clusters"
DimPlot(cntl, label = TRUE, reduction = "umap", pt.size = 1, repel = TRUE)
DimPlot(cntl, label = TRUE, reduction = "umap", pt.size = 1, repel = TRUE)


FeaturePlot(cntl, features = "ACE2", pt.size = 2, order = TRUE, split.by = "orig.ident") & theme(legend.position = "right") + theme_void()
FeaturePlot(cntl, features = "Ace2", pt.size = 2, order = TRUE, no.axes = TRUE, no.legend = FALSE, split.by = "orig.ident")

VlnPlot(cntl, features = "Ace2")
data <- data.frame(FetchData(object = cntl, slot = 'data',  vars = c('Ace2', 'ident', 'Condition'), group.by = "Condition"))
is.data.frame(data)
export(data, file = "expressionACE2.xlsx", rowNames = T, colNames = T)

FeaturePlot(cntl, features = "Cxcl16", pt.size = 2, order = TRUE) & theme(legend.position = "right") + theme_void()
FeaturePlot(cntl, features = "Cxcl16", pt.size = 2, order = TRUE, split.by = "orig.ident") & theme(legend.position = "right") + theme_void()
data <- data.frame(FetchData(object = cntl, slot = 'data',  vars = c('Cxcl16', 'ident', 'Condition'), group.by = "Condition"))
is.data.frame(data)
export(data, file = "expressionCxcl16 day0.xlsx", rowNames = T, colNames = T)

FeaturePlot(cntl, features = "Ccl17", pt.size = 2, order = TRUE, split.by = "orig.ident") & theme(legend.position = "right") + theme_void()
data <- data.frame(FetchData(object = cntl, slot = 'data',  vars = c('Ccl17', 'ident', 'Condition'), group.by = "Condition"))
is.data.frame(data)
export(data, file = "expressionCcl17 day0.xlsx", rowNames = T, colNames = T)

FeaturePlot(cntl, features = "Arg1", pt.size = 2, order = TRUE, split.by = "orig.ident")
FeaturePlot(cntl, features = "S100a4", pt.size = 0.7, order = TRUE, split.by = "orig.ident")
FeaturePlot(cntl, features = "Ifitm3", pt.size = 0.7, order = TRUE, split.by = "orig.ident")
FeaturePlot(cntl, features = "Apoe", pt.size = 0.7, order = TRUE, split.by = "orig.ident")
FeaturePlot(cntl, features = "C1qa", pt.size = 0.7, order = TRUE, split.by = "orig.ident")
FeaturePlot(cntl, features = "C1qb", pt.size = 0.7, order = TRUE, , split.by = "orig.ident")

#####Interstitialmacrophages
FeaturePlot(cntl, features = "Cxcl16", pt.size = 1, order = TRUE)
FeaturePlot(cntl, features = "C1qb", pt.size = 0.7, order = TRUE, split.by = "orig.ident")
FeaturePlot(cntl, features = "Apoe", pt.size = 0.7, order = TRUE, split.by = "orig.ident")
FeaturePlot(cntl, features = "C1qc", pt.size = 0.7, order = TRUE, split.by = "orig.ident")
FeaturePlot(cntl, features = "Cx3cr1", pt.size = 0.7, order = TRUE, split.by = "orig.ident")
FeaturePlot(cntl, features = "H2-Ab1", pt.size = 0.7, order = TRUE, , split.by = "orig.ident")
FeaturePlot(cntl, features = "Mafb", pt.size = 0.7, order = TRUE, , split.by = "orig.ident")
FeaturePlot(cntl, features = "Ccr2", pt.size = 0.7, order = TRUE, , split.by = "orig.ident")
FeaturePlot(cntl, features = "Apoe", pt.size = 0.7, order = TRUE, , split.by = "orig.ident")
FeaturePlot(cntl, features = "Mgl2", pt.size = 0.7, order = TRUE, , split.by = "orig.ident")

#####Alveolarmacrophages
FeaturePlot(cntl, features = "Marco", pt.size = 1, order = TRUE, split.by = "orig.ident")
FeaturePlot(cntl, features = "Mrc1", pt.size = 0.7, order = TRUE, split.by = "orig.ident")
FeaturePlot(cntl, features = "Chil3", pt.size = 0.7, order = TRUE, split.by = "orig.ident")
FeaturePlot(cntl, features = "Car4", pt.size = 0.7, order = TRUE, split.by = "orig.ident")
FeaturePlot(cntl, features = "Ear1", pt.size = 0.7, order = TRUE, split.by = "orig.ident")
FeaturePlot(cntl, features = "Ear2", pt.size = 0.7, order = TRUE, , split.by = "orig.ident")
FeaturePlot(cntl, features = "Plet1", pt.size = 0.7, order = TRUE, , split.by = "orig.ident")
FeaturePlot(cntl, features = "Fabp1", pt.size = 0.7, order = TRUE, , split.by = "orig.ident")
FeaturePlot(cntl, features = "Fabp4", pt.size = 0.7, order = TRUE, , split.by = "orig.ident")
FeaturePlot(cntl, features = c("Marco","Mrc1", "Chil3", "Car4"), pt.size = 0.7, order = TRUE)
FeaturePlot(cntl, features = c("Ear1", "Ear2", "Plet1", "Fabp1"), pt.size = 0.7, order = TRUE)
FeaturePlot(cntl, features = c("Ear1", "Ear2", "Plet1","Fabp4"), pt.size = 0.7, order = TRUE)
FeaturePlot(cntl, features = c("Marco","Mrc1", "Chil3", "Car4", "Ear1", "Ear2", "Plet1", "Fabp1", "Fabp4"), pt.size = 0.7, order = TRUE)

###Type 2 Cytokines and regulatory markers
FeaturePlot(cntl, features = "Il4", pt.size = 0.7, order = TRUE, , split.by = "orig.ident")
FeaturePlot(cntl, features = "Il5", pt.size = 0.7, order = TRUE, , split.by = "orig.ident")
FeaturePlot(cntl, features = "Il13", pt.size = 0.7, order = TRUE, , split.by = "orig.ident")
FeaturePlot(cntl, features = "Il2", pt.size = 0.7, order = TRUE, , split.by = "orig.ident")
FeaturePlot(cntl, features = "Il10", pt.size = 0.7, order = TRUE, , split.by = "orig.ident")
FeaturePlot(cntl, features = "Tgfb1", pt.size = 0.7, order = TRUE, , split.by = "orig.ident")
FeaturePlot(cntl, features = "Il27", pt.size = 0.7, order = TRUE, , split.by = "orig.ident")
FeaturePlot(cntl, features = "Foxp3", pt.size = 0.7, order = TRUE, , split.by = "orig.ident")


FeaturePlot(cntl, features = "Chil3", pt.size = 0.7, order = TRUE, , split.by = "orig.ident")
FeaturePlot(cntl, features = "Plet1", pt.size = 0.7, order = TRUE, , split.by = "orig.ident")
FeaturePlot(cntl, features = "Mrc1", pt.size = 0.7, order = TRUE, , split.by = "orig.ident")
FeaturePlot(cntl, features = "Siglecf", pt.size = 0.7, order = TRUE, , split.by = "orig.ident")
FeaturePlot(cntl, features = "Itgax", pt.size = 0.7, order = TRUE, , split.by = "orig.ident")

#####Tissue Repair
FeaturePlot(cntl, features = "Areg", pt.size = 1, order = TRUE, split.by = "orig.ident")
FeaturePlot(cntl, features = "Arg1", pt.size = 0.7, order = TRUE, split.by = "orig.ident")
FeaturePlot(cntl, features = "Chil3", pt.size = 0.7, order = TRUE, split.by = "orig.ident")
FeaturePlot(cntl, features = "Retnla", pt.size = 0.7, order = TRUE, split.by = "orig.ident")
FeaturePlot(cntl, features = "Car4", pt.size = 0.7, order = TRUE, , split.by = "orig.ident")

#Testing singleR main labels
library(magrittr)
library(ggplot2)
Idents(cntl) <- 'singleR.fine'
umap_fine = cntl@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% cbind(fine = cntl@meta.data$singleR.fine)


Idents(cntl) <- 'singleR.fine'
umap_fine = cntl@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% cbind(fine = cntl@meta.data$singleR.fine)

ggplot(umap_fine, aes(x=UMAP_1, y=UMAP_2, color=fine)) + geom_point() + 
  scale_color_manual(values=c("DC (DC.103-11B+24+)" = "blue", 
                              "DC (DC.103-11B+F4-80LO.KD)" = "blue", 
                              "DC (DC.103+11B-)" = "blue", 
                              "DC (DC.11B+)" = "blue", 
                              "DC (DC.8-)" = "blue", 
                              "DC (DC.8-4-11B-)" = "blue", 
                              "DC (DC.8-4-11B+)" = "blue",
                              "DC (DC.8+)" = "blue",
                              "DC (DC.LC)" = "blue",
                              "DC (DC.PDC.8-)" = "blue",
                              "DC (DC.PDC.8+)" = "blue",
                              "DC (DC)" = "blue",
                              "Epithelial cells (Ep.5wk.MEClo)" = "green2",
                              "Epithelial cells (Ep.8wk.CEC.Sca1+)" = "green2",
                              "Endothelial cells (BEC)" = "red",
                              "Endothelial cells (LEC.CFA)" = "red",
                              "Endothelial cells (LEC)" = "red",
                              "Stem cells (GMP)" = "yellow",
                              "Stem cells (LTHSC)" = "yellow",
                              "Stem cells (SC.CDP)" = "yellow",
                              "Stem cells (SC.LT34F)" = "yellow",
                              "Stem cells (SC.STSL)" = "yellow",
                              "Stem cells (SC.MEP)" = "yellow",
                              "Epithelial cells (MECHI.GFP-.ADULT)" = 
                                "green2",
                             "Epithelial cells (Ep.5wk.MEC.Sca1+)" = "green2",
                             "Epithelial cells (EP.MECHI)" = "green2",
                             "Epithelial cells (MECHI.GFP-.ADULT)" = "green2",
                             "Epithelial cells (MECHI.GFP+.ADULT.KO)" = "green2",
                             "Epithelial cells (MECHI.GFP+.ADULT.KO)" = "green2",
                             "Epithelial cells (MECHI.GFP+.ADULT)" = "green2"))

ggplot(umap_fine, aes(x=UMAP_1, y=UMAP_2, color=fine)) + geom_point() + 
  scale_color_manual(values=c("T cells (T.4FP3+25+)" = "blue", 
                              "T cells (T.4MEM44H62L)" = "blue",
                              "T cells (T.4MEM49D+11A+.D30.LCMV)" = "blue",
                              "T cells (T.4Nve)" = "blue",
                              "T cells (T.8EFF.OT1.12HR.LISOVA)" = "blue",
                              "T cells (T.8EFF.OT1.48HR.LISOVA)" = "blue",
                              "T cells (T.8EFF.OT1.D10LIS)" = "blue",
                              "T cells (T.8MEM.OT1.D106.VSVOVA)" = "blue",
                              "T cells (T.T.8MEM.OT1.D45.LISOVA)" = "blue",
                              "T cells (T.T.8Mem)" = "blue",
                              "T cells (T.8MEMKLRG1-CD127+.D8.LISOVA)" = "blue",
                              "T cells (T.8NVE.OT1)" = "blue",
                              "T cells (T.8Nve)" = "blue",
                              "T cells (T.8NVE)" = "blue",
                              "T cells (T.CD4.1H)" = "blue",
                              "T cells (T.CD4.24H)" = "blue",
                              "T cells (T.CD4.5H)" = "blue",
                              "T cells (T.CD4.CTR)" = "blue",
                              "T cells (T.CD4TESTCJ)" = "blue",
                              "T cells (T.CD8.1H)" = "blue",
                              "T cells (T.CD8.CTR)" = "blue",
                              "T cells (T.ETP)" = "blue",
                              "T cells (T.Tregs)" = "blue",
                              "ILC (ILC2)" = "orange",
                              "ILC (LIV.ILC1.DX5-)" = "orange",
                              "ILC (LIV.NK.DX5+)" = "orange",
                              "ILC (LPL.NCR+CNK)" = "orange",
                              "ILC (LPL.NCR+ILC1)" = "orange",
                               "NK cells (NK.49CI+)" = "red",
                              "NK cells (NK.49H+)" = "red",
                              "NK cells (NK.DAP10-)" = "red",
                              "NK cells (NK.H+MCMV1)" = "red",
                              "NK cells (NK.MCMV7)" = "red",
                              "NK cells (NK)" = "red",
                              "NKT (NKT.4-)" = "yellow",
                              "NKT (NKT.4+)" = "yellow",
                              "NKT (NKT.44+NK1.1-)" = "yellow",
                              "NKT (NKT.44+NK1.1+)" = "yellow",
                              "Tgd (Tgd.mat.VG1+VD6+)" = "green2",
                               "Tgd (Tgd.VG2+)" = "green2",
                              "Tgd (Tgd.mat.VG2+)" = "green2"))




ggplot(umap_fine, aes(x=UMAP_1, y=UMAP_2, color=fine)) + geom_point() + 
  scale_color_manual(values=c("Fibroblasts (FI.MTS15+)" = "blue", 
                              "Fibroblasts (FI)" = "blue",
                              "Fibroblasts (FRC)" = "blue",
                              "ILC (ILC2)" = "red",
                              "ILC (LIV.ILC1.DX5-)" = "red",
                              "ILC (LIV.NK.DX5+)" = "red",
                              "ILC (LPL.NCR+CNK)" = "red",
                              "ILC (LPL.NCR+ILC1)" = "red",
                              "Macrophages (MF.103-11B+24-)" = "green2",
                              "Macrophages (MF.11C-11B+)" = "green2",
                              "Macrophages (MF.480HI.NAIVE)" = "green2",
                              "Macrophages (MF.11CLOSER.SALM3)" = "green2",
                              "Macrophages (MF.480INT.NAIVE)" = "green2",
                              "Macrophages (MF.ALV)" = "red",
                              "Macrophages (MF.II-480HI)" = "green2",
                              "Macrophages (MF.II+480LO)" = "green2",
                              "Macrophages (MF)" = "green2",
                              "Macrophages (MF.MEDL)" = "green2",
                              "Macrophages (MF.SBCAPS)" = "green2",
"Macrophages (MF)" = "green2",
"Macrophages (MFIO5.II+480LO)" = "green2",
"Macrophages (MFIO5.II-480HI)" = "green2",
"Macrophages (MFIO5.II-480INT)" = "green2",
"Macrophages (MFIO5.II+480INT)" = "green2",
"Monocytes (MO.6C-II+)" = "yellow",
"Monocytes (MO.6C-IIINT)" = "yellow",
"Monocytes (MO.6C+II-)" = "yellow",
"Monocytes (MO.6C+II+)" = "yellow",
"Monocytes (MO)" = "yellow",
"Monocytes (MO.6C-II-)" = "yellow"))



ggplot(umap_fine, aes(x=UMAP_1, y=UMAP_2, color=fine)) + geom_point() + 
  scale_color_manual(values=c(
                              "Stromal cells (ST.31-38-44-)" = "blue",
                              "Stromal cells (DN.CFA)" = "blue",
                              "Stromal cells (DN)" = "blue",
                              "Neutrophils (GN.URAC)" = "red",
                              "Neutrophils (GN.Thio)" = "red",
                              "Neutrophils (GN.ARTH)" = "red",
                              "Monocytes (MO.6+2+)" = "green2",
                              "Monocytes (MO.6C-II+)" = "green2",
                              "Monocytes (MO.6C+II-)" = "green2",
                              "Monocytes (MO.6C+II+)" = "green2", 
                              "Fibroblasts (FI.MTS15+)" = "orange",
                              "Fibroblasts (FRC.CAD11.WT)" = "orange",
                              "Fibroblasts (FRC.CFA)" = "orange",
                              "Fibroblasts (FRC)" = "orange",
                              "Fibroblasts (FI)" = "orange"))



ggplot(umap_fine, aes(x=UMAP_1, y=UMAP_2, color=fine)) + geom_point() + 
  scale_color_manual(values=c(
    "Mast cells (MC.TR)" = "blue",
    "B cells (B.CD19CONTROL)" = "red",
    "B cells (B.Fo)" = "red",
    "B cells (B.FO)" = "red",
    "B cells (B.FrE)" = "red",
    "B cells (B.FRE)" = "red",
    "B cells (B.FrF)" = "red",
    "B cells (B.GC)" = "red",
    "B cells (B.MZ)" = "red",
    "B cells (B.T1)" = "red",
    "B cells (B.T2)" = "red",
    "B cells (B.T3)" = "red",
    "B cells (B1a)" = "red",
    "B cells (B1A)" = "red",
    "B cells (B1b)" = "red",
    "B cells (preB.FrD)" = "red",
    "B cells (proB.FRA)" = "red",
    "Basophils (BA)" = "green2",
"Eosinophils (EO)" = "yellow",
"Fibroblasts (FI.MTS15+)" = "orange",
"Fibroblasts (FRC.CAD11.WT)" = "orange",
"Fibroblasts (FRC.CFA)" = "orange",
"Fibroblasts (FRC)" = "orange",
"Fibroblasts (FI)" = "orange"))

ggplot(umap_fine, aes(x=UMAP_1, y=UMAP_2, color=fine)) + geom_point() + 
  scale_color_manual(values=c(
    "NK cells (NK.DAP10-)" = "blue",
    "B cells (B.CD19CONTROL)" = "red",
    "B cells (B.Fo)" = "red",
    "B cells (B.FO)" = "red",
    "B cells (B.FrE)" = "red",
    "B cells (B.FRE)" = "red",
    "B cells (B.FrF)" = "red",
    "B cells (B.GC)" = "red",
    "B cells (B.MZ)" = "red",
    "B cells (B.T1)" = "red",
    "B cells (B.T2)" = "red",
    "B cells (B.T3)" = "red",
    "B cells (B1a)" = "red",
    "B cells (B1A)" = "red",
    "B cells (B1b)" = "red",
    "B cells (preB.FrD)" = "red",
    "B cells (proB.FRA)" = "red",
    "Basophils (BA)" = "green2",
    "Eosinophils (EO)" = "yellow",
    "Neutrophils (GN.ARTH)" = "orange",
    "Neutrophils (GN.Thio)" = "orange",
    "Neutrophils (GN.URAC)" = "orange",
    "Neutrophils (GN)" = "orange"))






#Rename cell clusters
DefaultAssay(cntl) <- "RNA"
Idents(cntl) <- "seurat_clusters"
current.cluster.ids <- c (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19)
new.cluster.ids <- c("Endothelial cells", "Endothelial cells", "Fibroblast-stromal cells", "B cells", "Fibroblast-stromal cells", "Other Lymphocytes", "Fibroblast-stromal cells", "Dendritic cells", 
                     "Mono-Macs", "Endothelial cells", "Alv.Macs", "AA Alv.Macs", "Endothelial cells", "Endothelial cells", "Endothelial cells", "Granulocytes", "Epithelial cells", "NK cells", "Fibroblast-stromal cells", "Brc5+ Alv.Macs")

names(new.cluster.ids) <- levels(cntl)
cntl <- RenameIdents(cntl, new.cluster.ids)
DimPlot(cntl, label = TRUE, reduction = "umap", pt.size = 1, repel = TRUE)
DimPlot(cntl, label = TRUE, reduction = "umap", pt.size = 0.7, repel = TRUE) +  theme_void() + NoLegend() 
ggsave("Seuratclusters.svg", width = 10, height = 10, dpi=700)
cntl[["CellType"]] <- Idents(object = cntl)


# Rename identity classes
Idents(cntl) <- "orig.ident"
new.cluster.id <- c("nipo", "control")
names(new.cluster.id) <- levels(cntl)
cntl <- RenameIdents(cntl, "control" = "control", "nipo" = "nipo")
DimPlot(cntl, label = TRUE, reduction = "umap", pt.size = 1, repel = TRUE)
cntl[["Condition"]] <- Idents(object = cntl)

Idents(cntl) <- "seurat_clusters"
Idents(cntl) <- "CellType"

Idents(cntl) <- "CellType"
new.idents <- cntl@meta.data$CellType #add cell type as idents
Idents(object = cntl) <- new.idents #add cell type as idents





DimPlot(cntl, label = TRUE, reduction = "umap", pt.size = 1, repel = TRUE, split.by = "Condition")
DimPlot(cntl, label = FALSE, reduction = "umap", pt.size = 1, repel = TRUE) + NoLegend () + theme_void()
ggsave("Seuratclusters_Condition.svg", width = 10, height = 10, dpi=700)

cntl.markers <- FindAllMarkers(cntl, only.pos = T)
cntl.markers <- cntl.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
DoHeatmap(cntl, cntl.markers$gene, size = 5, label = FALSE)+
  scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#6300ae','#F8766D','#FFCC99')), midpoint = 0, guide = "colourbar", aesthetics = "fill") 
cntl.rna.markers <- cntl.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
export(cntl.rna.markers, file = "CORRECTED_RNAmarkers_top50.xlsx")


Idents(cntl) <- "orig.ident"

saveRDS(cntl, file = "cntlclusters.rds")
cntl <- readRDS(file = "cntlclusters.rds")


Idents(cntl) <- "CellType"
# How many cells are in each cluster
table(Idents(cntl))

# How many cells are in each replicate?
table(cntl$Condition)

# What proportion of cells are in each cluster?
prop.table(table(Idents(cntl)))

# How does cluster membership vary by replicate?
table(Idents(cntl), cntl$Condition)

prop.table <- prop.table(table(Idents(cntl), cntl$Condition), margin = 2)

prop.table <-data.frame((prop.table))

dittoBarPlot(
  object = cntl,
  var = "CellType",
  group.by = "Condition",
  do.hover = T,
  ####color.panel = c("red", "orange", "purple"),
  theme =theme_classic( base_size = 15
  ))

ggsave("DittoProportion of cells.svg", width = 10, height = 10, dpi=700)

dittoBarPlot(
  object = cntl,
  var = "CellType",
  group.by = "Condition",
  do.hover = T,
  theme =theme_classic( base_size = 15
  ) + NoLegend()) 

#AddModuleScore can be used to see if expression of given gene set is 
#enriched vs set of randomly selected (but based on expression bins) control genes. 
#This might help to clean up the plot as it sounds like the enrichment of the whole gene set 
#would likely be cell type specific whereas one particular gene 
#might also be expressed in other cell types.

DefaultAssay(cntl) <- "RNA"
Idents(cntl) <- "singleR.main"
Idents(cntl) <- "CellType"
Idents(cntl) <- "orig.ident"
Idents(cntl) <- "Condition"
Type2_cytokine_gene_list <- list(c("Il4", "Il13"))
cntl <- AddModuleScore(cntl, features = Type2_cytokine_gene_list, name = "Type2genes_score")
names(x = cntl[[]])
FeaturePlot(object = cntl, features = "Type2genes_score1", pt.size = 1, order = TRUE, split.by = "orig.ident")
DotPlot(object = cntl, features = "Type2genes_score1", split.by = "orig.ident")

cytokines <- c("Il4", "Il13", "Il5", "Il10")
###playing around with Dot plots
DotPlot(object = cntl, features =  "Type2genes_score1", cols = c("lightgrey", "blue"), split.by = "orig.ident")
DotPlot(object = cntl, features = cytokines, scale = TRUE)
ggsave("Type 2 cytokinesclusters.svg", width = 10, height = 10, dpi=700)
# limits expression to only a cell type of interest
DotPlot(object = cntl, features = cytokines, idents = "Granulocytes", split.by = "orig.ident", scale = TRUE) 
VlnPlot(object = cntl, features = "Type2genes_score1", split.by = "orig.ident")


#Differential Gene Expression
DefaultAssay(cntl) <- "RNA"
Idents(cntl) <- "singleR.main"
Idents(cntl) <- "CellType"
Idents(cntl) <- "CellType"
Idents(cntl) <- "orig.ident"
Idents(cntl) <- "Condition"
Idents(cntl) <- "seuratclusters"

Idents(cntl) <- "CellType"
Alv.Macs <- subset(cntl, idents = "Alv.Macs")
saveRDS(Alv.Macs, file = "cntlAlv.Macs.rds")
Idents(Alv.Macs) <- "orig.ident"
avg.Alv.Mcas <- log1p(AverageExpression(Alv.Macs, verbose = FALSE)$RNA)
head(avg.Alv.Mcas, n = 20)
avg.Alv.Mcas <-data.frame(avg.Alv.Mcas)
is.data.frame(avg.Alv.Mcas)
export(avg.Alv.Mcas, file = "avg.Alv.Mcas.xlsx", rowNames = T, colNames = T)
Alv.Macs.de.markers <- FindMarkers(Alv.Macs, ident.1 = "nipo", ident.2 = "control")
head(Alv.Macs.de.markers, n = 300)
export(Alv.Macs.de.markers, file = "Alv.Macs.xlsx", rowNames = T, colNames = T)

Idents(cntl) <- "CellType"
Version(cntl)
Alv.Macs <- subset(cntl, idents = "Alv.Macs")
saveRDS(Alv.Macs, file = "cntlAlv.macs.rds")
Idents(Alv.Macs) <- "orig.ident"
avg.Alv.Macs <- log1p(AverageExpression(Alv.Macs, verbose = FALSE)$RNA)
head(avg.Alv.Macs, n = 20)
avg.Alv.Macs <-data.frame(avg.Alv.Macs)
is.data.frame(avg.Alv.Macs)
export(avg.Alv.Macs, file = "avg.Alv.Macs.xlsx", rowNames = T, colNames = T)
Alv.Macs.de.markers <- FindMarkers(Alv.Macs, ident.1 = "nipo", ident.2 = "control")
head(Alv.Macs.de.markers, n = 300)
export(Alv.Macs.de.markers, file = "Alv.Macsmarkersnippovscontrol.xlsx", rowNames = T, colNames = T)

Idents(cntl) <- "CellType"
cntl.myeloid <- subset(cntl, idents = c("Alv.Macs", "AA Alv.Macs", "Mono-Macs", "Brc5+ Alv.Macs", "Dendritic cells"))
cntl.mono.macs <- subset(cntl, idents = c("Alv.Macs", "AA Alv.Macs", "Mono-Macs", "Brc5+ Alv.Macs"))
saveRDS(cntl.myeloid, file = "cntl.myeloid.rds")
saveRDS(cntl.mono.macs, file = "cntl.mono.macs.rds")
Idents(cntl.myeloid) <- "orig.ident"
avg.cntl.myeloid <- log1p(AverageExpression(cntl.myeloid, verbose = FALSE)$RNA)
head(avg.cntl.myeloid, n = 50)
avg.cntl.myeloid <-data.frame(avg.cntl.myeloid)
is.data.frame(avg.cntl.myeloid)
export(avg.cntl.myeloid, file = "avg.cntl.myeloid.xlsx", rowNames = T, colNames = T)
Myeloid.de.markers <- FindMarkers(cntl.myeloid, ident.1 = "nipo", ident.2 = "control")
head(Myeloid.de.markers, n = 300)
export(Myeloid.de.markers, file = "Myeloidmarkersnippovscontrol.xlsx", rowNames = T, colNames = T)

Idents(cntl) <- "CellType"
cntl.AMs <- subset(cntl, idents = c("Alv.Macs", "AA Alv.Macs"))
saveRDS(cntl.AMs, file = "cntl.AMs2023.rds")
Idents(cntl.AMs) <- "orig.ident"
avg.cntl.AMs <- log1p(AverageExpression(cntl.AMs, verbose = FALSE)$RNA)
head(avg.cntl.AMs, n = 50)
avg.cntl.AMs <-data.frame(avg.cntl.AMs)
is.data.frame(avg.cntl.AMs)
export(avg.cntl.AMs, file = "avg.cntl.AMs.xlsx", rowNames = T, colNames = T)
AMs.de.markers <- FindMarkers(cntl.AMs, ident.1 = "nipo", ident.2 = "control")
head(AMs.de.markers, n = 300)
export(AMs.de.markers, file = "AMmarkersnippovscontrol.xlsx", rowNames = T, colNames = T)


Idents(cntl) <- "CellType"
cntl.AMs <- subset(cntl, idents = c("Alv.Macs", "AA Alv.Macs", "Brc5+ Alv.Macs"))
saveRDS(cntl.AMs, file = "cntl.AMs+Brc52023.rds")
Idents(cntl.AMs) <- "orig.ident"
avg.cntl.AMs <- log1p(AverageExpression(cntl.AMs, verbose = FALSE)$RNA)
head(avg.cntl.AMs, n = 50)
avg.cntl.AMs <-data.frame(avg.cntl.AMs)
is.data.frame(avg.cntl.AMs)
export(avg.cntl.AMs, file = "avg.cntl.AMs.xlsx", rowNames = T, colNames = T)
AMs.de.markers <- FindMarkers(cntl.AMs, ident.1 = "nipo", ident.2 = "control")
head(AMs.de.markers, n = 300)
export(AMs.de.markers, file = "AMmarkersnippovscontrol.xlsx", rowNames = T, colNames = T)


Idents(cntl) <- "CellType"
AAAlv.Macs <- subset(cntl, idents = "AA Alv.Macs")
saveRDS(AAAlv.Macs, file = "cntlAAAlv.macs.rds")
Idents(AAAlv.Macs) <- "orig.ident"
avg.AAAlv.Macs <- log1p(AverageExpression(AAAlv.Macs, verbose = FALSE)$RNA)
head(avg.AAAlv.Macs, n = 40)
avg.AAAlv.Macs <-data.frame(avg.AAAlv.Macs)
is.data.frame(avg.AAAlv.Macs)
export(avg.AAAlv.Macs, file = "avg.AAAlv.Macs.xlsx", rowNames = T, colNames = T)
AAAlv.Macs.de.markers <- FindMarkers(AAAlv.Macs, ident.1 = "nipo", ident.2 = "control")
head(AAAlv.Macs.de.markers, n = 300)
export(AAAlv.Macs.de.markers, file = "AAAlv.Macsmarkersnippovscontrol.xlsx", rowNames = T, colNames = T)

Idents(cntl) <- "CellType"
AAAlv.Macs <- subset(cntl, idents = "AA.AlvMacs")
saveRDS(AAAlv.Macs, file = "cntlAAAlv.macs.rds")
Idents(AAAlv.Macs) <- "orig.ident"
avg.AAAlv.Macs <- log1p(AverageExpression(AAAlv.Macs, verbose = FALSE)$RNA)
head(avg.AAAlv.Macs, n = 40)
avg.AAAlv.Macs <-data.frame(avg.AAAlv.Macs)
is.data.frame(avg.AAAlv.Macs)
export(avg.AAAlv.Macs, file = "avg.AAAlv.Macs.xlsx", rowNames = T, colNames = T)
AAAlv.Macs.de.markers <- FindMarkers(AAAlv.Macs, ident.1 = "nipo", ident.2 = "control")
head(AAAlv.Macs.de.markers, n = 300)
export(AAAlv.Macs.de.markers, file = "AAAlv.Macsmarkersnippovscontrol.xlsx", rowNames = T, colNames = T)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")

browseVignettes("EnhancedVolcano")

####Interpretation of results
####p_val : p_val (unadjusted)
###avg_logFC : log fold-chage of the average expression between the two groups. Positive values indicate that the feature is more highly expressed in the first group.
###pct.1 : The percentage of cells where the feature is detected in the first group
####pct.2 : The percentage of cells where the feature is detected in the second group
#####p_val_adj : Adjusted p-value, based on bonferroni correction using all features in the dataset.

# Pre-filter features that are detected at <50% frequency in either CD14+ Monocytes or FCGR3A+
# Monocytes
###head(FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono", min.pct = 0.5))

# Pre-filter features that have less than a two-fold change between the average expression of
# CD14+ Monocytes vs FCGR3A+ Monocytes
###head(FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono", logfc.threshold = log(2)))



Idents(cntl) <- "CellType"
Tcells <- subset(cntl, idents = "Other Lymphocytes")
saveRDS(Tcells, file = "cntlotherlymphocytes.rds")
Idents(Tcells) <- "orig.ident"
avg.Tcells <- log1p(AverageExpression(Tcells, verbose = FALSE)$RNA)
head(avg.Tcells, n = 20)
avg.Tcells <-data.frame(avg.Tcells)
is.data.frame(avg.Tcells)
export(avg.Tcells, file = "Avg.T cells.xlsx", , rowNames = T, colNames = T)
Tcells.de.markers <- FindMarkers(Tcells, ident.1 = "nipo", ident.2 = "control")
head(Tcells.de.markers, n = 300)
export(Tcells.de.markers, file = "Tcellssmarkersnippovscontrol.xlsx", rowNames = T, colNames = T)


#Plot2
genes.to.label = c("Ly6a", "Ly6c1", "Fabp5", "Alox5ap", "Ms4a4c", "Cldn10", "Gbp2", "Gbp7", "Pdcd1lg2", "Jun")
p5 <- ggplot(avg.Macrophages, aes(nipo, control)) + geom_point(color = "steelblue2") + ggtitle("Macrophages")
p5 <- LabelPoints(plot = p5, points = genes.to.label, repel = TRUE, max.overlaps = Inf)
###p6 <- ggplot(avg.Mono.Mac, aes(nippocov2, controlcov2)) + geom_point(color = "steelblue2") + ggtitle("Mono.Mac")
###p6 <- LabelPoints(plot = p6, points = genes.to.label, repel = TRUE, max.overlaps = Inf)
p7 <- ggplot(avg.Tcells, aes(nippocov2, controlcov2)) + geom_point(color = "steelblue2") + ggtitle("T")
p7 <- LabelPoints(plot = p7, points = genes.to.label, repel = TRUE, max.overlaps = Inf)
p8 <- ggplot(avg.Stromal.Fibroblasts, aes(nippocov2, controlcov2)) + geom_point(color = "steelblue2") + ggtitle("Stromal.Fibroblasts")
plot_grid(p5)

Idents(cntl) <- "CellType"
Monocytes <- subset(cntl, idents = "Mono-Macs")
saveRDS(Tcells, file = "cntlMonocytes.rds")
Idents(Monocytes) <- "orig.ident"
avg.Monocytes <- log1p(AverageExpression(Monocytes, verbose = FALSE)$RNA)
head(avg.Monocytes, n = 20)
avg.Monocytes <-data.frame(avg.Monocytes)
is.data.frame(avg.Monocytes)
export(avg.Monocytes, file = "Avg.Monocytes.xlsx")

Idents(cntl) <- "CellType"
Neutrophils <- subset(cntl, idents = "Neutrophils")
Idents(Neutrophils) <- "orig.ident"
avg.Neutrophils <- log1p(AverageExpression(Neutrophils, verbose = FALSE)$RNA)
head(avg.Neutrophils, n = 20)
avg.Neutrophils <-data.frame(avg.Neutrophils)
is.data.frame(avg.Neutrophils)
export(avg.Neutrophils, file = "Avg.Neutrophils.xlsx")


Endothelial <- subset(cntl, idents = "Endothelial cells")
Idents(Endothelial) <- "orig.ident"
avg.Endothelial <- log1p(AverageExpression(Endothelial, verbose = FALSE)$RNA)
head(avg.Endothelial, n = 20)
avg.Endothelial <-data.frame(avg.Endothelial)
is.data.frame(avg.Endothelial)
export(avg.Endothelial, file = "Avg.Endothelial.xlsx")


Idents(cntl) <- "singleR.main"
Stromal <- subset(cntl, idents = "Stromal cells")
Idents(Stromal) <- "orig.ident"
avg.Stromal <- log1p(AverageExpression(Stromal, verbose = FALSE)$RNA)
head(avg.Stromal, n = 20)
avg.Stromal <-data.frame(avg.Stromal)
is.data.frame(avg.Stromal)
export(avg.Stromal, file = "Avg.Stromal.xlsx")

#Plot2
genes.to.label = c("Ly6a", "Ly6c1", "Fabp5", "Alox5ap", "Ms4a4c", "Cldn10", "Gbp2", "Gbp7", "Pdcd1lg2", "Jun")
p5 <- ggplot(avg.Monocytes, aes(nipo, control)) + geom_point(color = "steelblue2") + ggtitle("Monocytes")
p5 <- LabelPoints(plot = p5, points = genes.to.label, repel = TRUE, max.overlaps = Inf)
p6 <- ggplot(avg.Neutrophils, aes(nipo, control)) + geom_point(color = "steelblue2") + ggtitle("Neutrophils")
p6 <- LabelPoints(plot = p6, points = genes.to.label, repel = TRUE, max.overlaps = Inf)
p7 <- ggplot(avg.Endothelial, aes(nipo, control)) + geom_point(color = "steelblue2") + ggtitle("Endothelial.1")
p7 <- LabelPoints(plot = p7, points = genes.to.label, repel = TRUE, max.overlaps = Inf)
p8 <- ggplot(avg.Stromal, aes(nipo, control)) + geom_point(color = "steelblue2") + ggtitle("Endothelial.2")
p8 <- LabelPoints(plot = p7, points = genes.to.label, repel = TRUE, max.overlaps = Inf)
plot_grid(p5, p6, p7, p8)

Idents(cntl) <- "singleR.main"
Epithelial <- subset(cntl, idents = "Epithelial cells")
Idents(Epithelial) <- "orig.ident"
avg.Epithelial <- log1p(AverageExpression(Epithelial, verbose = FALSE)$RNA)
head(avg.Epithelial, n = 20)
avg.Epithelial <-data.frame(avg.Epithelial)
is.data.frame(avg.Epithelial)
export(avg.Epithelial, file = "Avg.Epithelial.xlsx")



Idents(cntl) <- "singleR.main"
B <- subset(cntl, idents = "B cells")
Idents(B) <- "orig.ident"
avg.B <- log1p(AverageExpression(B, verbose = FALSE)$RNA)
head(avg.B, n = 20)
avg.B <-data.frame(avg.B)
is.data.frame(avg.B)
export(avg.B, file = "Avg.B.xlsx")

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


Idents(cntl) <- "CellType"
Transitional.Mac <- subset(cov2, idents = "Transitional.Mac")
Idents(Transitional.Mac) <- "Condition"
avg.Transitional.Mac <- log1p(AverageExpression(Transitional.Mac, verbose = FALSE)$RNA)
head(avg.Transitional.Mac, n = 20)
avg.Transitional.Mac <-data.frame(avg.Transitional.Mac)
is.data.frame(avg.Transitional.Mac)

Idents(cntl) <- "single.Rmain"
DC <- subset(cntl, idents = "DC")
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


Idents(cntl) <- "singleR.main"
C0C1dge <- FindMarkers(cntl, ident.1 = "Macrophages", ident.2 = "Monocytes", ident.3 = "DC", verbose = FALSE)
head(C0C1dge, n = 50)
table <- head(C0C1dge, n = 50)
table <-data.frame(table)
is.data.frame(table)
export(table, file = "C01C1dge Top 50 Genes.xlsx")
rownames(head(C0C1dge, n = 50))

#DGE of chow vs hp in Cluster 1, 2, 6, 9, 10, 14, 15, 16
Idents(cntl) <- "singleR.main"
Idents(cntl) <- "seurat_clusters"
Idents(Macrophages) <-"orig.ident"
zero <- subset(cntl, idents = "Macrophages")
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


library(dittoSeq)

dittoBarPlot(
  object = cntl,
  var = "CellType",
  group.by = "Condition",
  do.hover = T,
  theme =theme_classic( base_size = 15
  ))

