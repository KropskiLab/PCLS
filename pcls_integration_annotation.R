library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
options(future.globals.maxSize = 40000 * 1024^2)
library(RColorBrewer)
library(viridis)
library(hrbrthemes)
library(SeuratDisk)

#read in pcls data
pcls.data <- Read10X_h5("~/Dropbox/Lab data/PCLS/filtered_feature_bc_matrix.h5")
pcls <- CreateSeuratObject(counts = pcls.data, min.features = 750, project = "PCLS")

#read in fresh counts from same donor
hd105.data <- Read10X_h5("~/Downloads/F02609_filtered_feature_bc_matrix.h5")
hd105 <- CreateSeuratObject(counts = hd105.data, min.features = 750, project = "VU_HD_105")

#read in controls from Habermann Science Advances 2020
hd65.data <- Read10X_h5("~/Desktop/sa_control/vuhd065_filtered_feature_bc_matrix.h5")
hd66.data <- Read10X_h5("~/Desktop/sa_control/vuhd066_filtered_feature_bc_matrix.h5")
hd67.data <- Read10X_h5("~/Desktop/sa_control/vuhd067_filtered_feature_bc_matrix.h5")
hd68.data <- Read10X_h5("~/Desktop/sa_control/vuhd068_filtered_feature_bc_matrix.h5")
hd69.data <- Read10X_h5("~/Desktop/sa_control/vuhd069_filtered_feature_bc_matrix.h5")
hd70.data <- Read10X_h5("~/Desktop/sa_control/vuhd070_filtered_feature_bc_matrix.h5")
hd71.data <- Read10X_h5("~/Desktop/sa_control/vuhd071_filtered_feature_bc_matrix.h5")
thd001.data <- Read10X_h5("~/Desktop/sa_control/thd001_filtered_feature_bc_matrix.h5")
thd002.data <- Read10X_h5("~/Desktop/sa_control/thd002_filtered_feature_bc_matrix.h5")
thd005a.data <- Read10X_h5("~/Desktop/sa_control/thd005a_filtered_feature_bc_matrix.h5")
thd005b.data <- Read10X_h5("~/Desktop/sa_control/thd005b_filtered_feature_bc_matrix.h5")
thd005c.data <- Read10X_h5("~/Desktop/sa_control/thd005c_filtered_feature_bc_matrix.h5")

hd65 <- CreateSeuratObject(counts = hd65.data, min.features = 750, project = "VU_HD_65")
hd66 <- CreateSeuratObject(counts = hd66.data, min.features = 750, project = "VU_HD_66")
hd67 <- CreateSeuratObject(counts = hd67.data, min.features = 750, project = "VU_HD_67")
hd68 <- CreateSeuratObject(counts = hd68.data, min.features = 750, project = "VU_HD_68")
hd69 <- CreateSeuratObject(counts = hd69.data, min.features = 750, project = "VU_HD_69")
hd70 <- CreateSeuratObject(counts = hd70.data, min.features = 750, project = "VU_HD_70")
hd71 <- CreateSeuratObject(counts = hd71.data, min.features = 750, project = "VU_HD_71")
thd001 <- CreateSeuratObject(counts = thd001.data, min.features = 750, project = "T_HD_001")
thd002 <- CreateSeuratObject(counts = thd002.data, min.features = 750, project = "T_HD_002")
thd005a <- CreateSeuratObject(counts = thd005a.data, min.features = 750, project = "T_HD_005a")
thd005b <- CreateSeuratObject(counts = thd005b.data, min.features = 750, project = "T_HD_005b")
thd005c <- CreateSeuratObject(counts = thd005c.data, min.features = 750, project = "T_HD_005c")

#merge
combined <- merge(hd65, y=c(hd66, hd67, hd68, hd69, hd70, hd71, thd001, thd002, thd005a, thd005b, thd005c, hd105, pcls))

#add metadata
Idents(combined) <- 'orig.ident'
Idents(combined, cells = WhichCells(combined, idents = c('T_HD_001', 'T_HD_002', 'T_HD_005a', 'T_HD_005b', 'T_HD_005c', 'VU_HD_65', 'VU_HD_66', 'VU_HD_67', 'VU_HD_68', 'VU_HD_69', 'VU_HD_70', 'VU_HD_71', 'VU_HD_105'))) <- "Fresh Control"
Idents(combined, cells = WhichCells(combined, idents = c('PCLS'))) <- "PCLS"
combined$diagnosis <- Idents(combined)

combined <- PercentageFeatureSet(combined, pattern = "^MT-", col.name = "percent.mt")
combined <- subset(combined, subset = percent.mt < 15 & percent.mt >0.5 & nFeature_RNA >200)



#perform integration
combined.list <- SplitObject(combined, split.by = "diagnosis")
for (i in names(combined.list)) {
  combined.list[[i]] <- SCTransform(combined.list[[i]], verbose = FALSE)
}

combined.features <- SelectIntegrationFeatures(object.list = combined.list, nfeatures = 3000)

combined.list <- PrepSCTIntegration(object.list = combined.list, anchor.features = combined.features)

combined.anchors <- FindIntegrationAnchors(object.list = combined.list, normalization.method = "SCT", 
                                           anchor.features = combined.features)
combined.integrated <- IntegrateData(anchorset = combined.anchors, normalization.method = "SCT")

combined.integrated <- RunPCA(object = combined.integrated, verbose = FALSE)
combined.integrated <- RunUMAP(object = combined.integrated, dims = 1:35)

DimPlot(combined.integrated, group.by=c("diagnosis"))

DefaultAssay(combined.integrated) <- 'integrated'
combined.integrated <- FindNeighbors(combined.integrated, dims = 1:35)
combined.integrated <- FindClusters(combined.integrated, resolution = 0.8)
DimPlot(combined.integrated, label = T)

DefaultAssay(combined.integrated) <- 'SCT'

DotPlot(combined.integrated, features = c("PTPRC", "EPCAM", "PECAM1", "CD3E", "COL1A1", "CD19", 'ACTA2', 'WT1', 'MKI67'))+ theme(axis.text.x = element_text(angle = 45, hjust=1))

epi_markers <- c("PTPRC", "PECAM1", "COL1A1", "MKI67", "TP63", 'KRT5', "SCGB1A1", "SCGB3A1", "SCGB3A2", "MUC5B", "MUC5AC", "FOXJ1", "TP73", "LAMP3", "ABCA3", "CEACAM6", "AGER", 'HOPX', 'NKX2-1', 'CLDN18', 'CLDN4', "KRT17", 'KRT8', 'KRT19', 'CDKN1A', 'SFN', 'CDKN2A', 'WWTR1', 'COL4A4', 'FBLN5', 'ELN', 'FOXI1', 'CALCA', 'AXIN2', 'NTM', 'NCKAP5',  'RTKN2', 'SCEL', 'LAMA3','COL4A1', 'COL4A2', 'COL4A3', 'COL4A5', 'COL4A6')
stromal_genes <- c("PTPRC", "EPCAM", "PECAM1", 'HEY1', 'CA4', 'PLVAP', 'COL15A1', "ACKR1", 'VWF', 'APLN', 'APLNR', 'CCL21','ACTA2', 'CSPG4', 'COL1A1', 'COL3A1', 'DCN', 'PLIN2', 'HAS1', 'WNT5A', 'WNT2', 'TCF21', 'CTHRC1', 'PI16', 'PDGFRA', 'PDGFRB', 'WT1', "MKI67")
immune_genes <- c("PECAM1", "EPCAM", "COL1A1", "PTPRC", "LYZ", "CD68", 'PPARG', "CD14", "FCGR3A", "CD86", 'HLA-DRA', "SPP1", 'MRC1', 'IL1B', "CD3E", "CD3G", 'CD3D', "CD4", "CD8A", "PDCD1", "IL7R", "FOXP3", "GNLY", "NKG7", "CPA3", "CD19", "JCHAIN", "IRF7", "MKI67")

DotPlot(combined.integrated, features = epi_markers)+ theme(axis.text.x = element_text(angle = 45, hjust=1))
DotPlot(combined.integrated, features = stromal_genes)+ theme(axis.text.x = element_text(angle = 45, hjust=1))
DotPlot(combined.integrated, features = immune_genes)+ theme(axis.text.x = element_text(angle = 45, hjust=1))

#subcluster on cluster 24
FindSubCluster(combined.integrated, "24", "integrated_snn", subcluster.name = "sub.cluster",
               resolution = 0.5,
               algorithm = 1
)

#remove doublet clusters
DefaultAssay(combined.integrated) <- 'integrated'
combined.integrated <- subset(combined.integrated, idents = c(0:11,13:15,18:21))
combined.integrated <- RunPCA(object = combined.integrated, verbose = FALSE)
combined.integrated <- RunUMAP(object = combined.integrated, dims = 1:30)
combined.integrated <- FindNeighbors(combined.integrated, dims = 0:30)
combined.integrated <- FindClusters(combined.integrated, resolution = 1.2)
DimPlot(combined.integrated, group.by=c("Status"))
DimPlot(combined.integrated, label = T)
DotPlot(combined.integrated, features = c("PTPRC", "EPCAM", "PECAM1", "CD3E", "COL1A1", "CD19", 'ACTA2', 'WT1', 'MKI67'))

DefaultAssay(combined.integrated) <- 'SCT'
DotPlot(combined.integrated, features = c("PTPRC", "PECAM1", "COL1A1",  "EPCAM", 'SCGB1A1', 'SCGB3A1', "SCGB3A2", 'MUC5B', 'KRT5', 'TP63', "KRT17", "ABCA3", 'LAMP3', "AGER", 'COL4A3','CDKN1A', 'KRT8', 'CLDN4', "FOXJ1", "MKI67", "CALCA", "FOXI1"))+ theme(axis.text.x = element_text(angle = 45, hjust=1))

epi_markers <- c("EPCAM", "PTPRC", "PECAM1", "COL1A1", "MKI67", "TP63", 'KRT5', "SCGB1A1", "SCGB3A1", "SCGB3A2", "MUC5B", "MUC5AC", "FOXJ1", "TP73", "LAMP3", "ABCA3", "CEACAM6", "AGER", 'HOPX', 'NKX2-1', 'CLDN18', 'CLDN4', "KRT17", 'KRT8', 'KRT19', 'CDKN1A', 'SFN', 'CDKN2A', 'WWTR1', 'COL4A4', 'FBLN5', 'ELN', 'FOXI1', 'CALCA', 'AXIN2', 'NTM', 'NCKAP5',  'RTKN2', 'SCEL', 'LAMA3','COL4A1', 'COL4A2', 'COL4A3', 'COL4A5', 'COL4A6')
stromal_genes <- c("PTPRC", "EPCAM", "PECAM1", 'HEY1', 'CA4', 'PLVAP', 'COL15A1', "ACKR1", 'VWF', 'APLN', 'APLNR', 'CCL21','ACTA2', 'CSPG4', 'COL1A1', 'COL3A1', 'DCN', 'PLIN2', 'HAS1', 'WNT5A', 'WNT2', 'TCF21', 'CTHRC1', 'PI16', 'PDGFRA', 'PDGFRB', 'WT1', "MKI67")
immune_genes <- c("PECAM1", "EPCAM", "COL1A1", "PTPRC", "LYZ", "CD68", 'PPARG', "CD14", "FCGR3A", "CD86", 'HLA-DRA', "SPP1", 'MRC1', 'IL1B', "CD3E", "CD3G", 'CD3D', "CD4", "CD8A", "PDCD1", "IL7R", "FOXP3", "GNLY", "NKG7", "CPA3", "CD19", "JCHAIN", "IRF7", "MKI67")

DefaultAssay(adata) <- 'SCT'

DotPlot(adata, features = epi_markers)+ theme(axis.text.x = element_text(angle = 45, hjust=1))
DotPlot(adata, features = stromal_genes)+ theme(axis.text.x = element_text(angle = 45, hjust=1))
DotPlot(adata, features = immune_genes)+ theme(axis.text.x = element_text(angle = 45, hjust=1))

#annotate epithelial cells
DimPlot(adata,label = T)
Idents(adata) <- 'sub.cluster'
Idents(adata, cells = WhichCells(adata, idents = c("24_4"))) <- "Transitional"
Idents(adata, cells = WhichCells(adata, idents = c("24_2"))) <- "Secretory - MUC5AC+"
Idents(adata, cells = WhichCells(adata, idents = c("24_0"))) <- "Secretory - SCGB1A1+/MUC5B+"
Idents(adata, cells = WhichCells(adata, idents = c("24_3"))) <- "Secretory - SCGB1A1+/SCGB3A1+"
Idents(adata, cells = WhichCells(adata, idents = c("24_1"))) <- "Secretory - SCGB1A1+/SCGB3A2+"
Idents(adata, cells = WhichCells(adata, idents = c("4"))) <- "Ciliated"
Idents(adata, cells = WhichCells(adata, idents = c("24_5"))) <- "Basal"
Idents(adata, cells = WhichCells(adata, idents = c("9","12", "23", "25"))) <- "AT2"
Idents(adata, cells = WhichCells(adata, idents = c("13"))) <- "AT1"

#annotate stromal cells
Idents(adata, cells = WhichCells(adata, idents = c("32"))) <- "Pericyte/SMC"
Idents(adata, cells = WhichCells(adata, idents = c("26"))) <- "FB"
Idents(adata, cells = WhichCells(adata, idents = c("28"))) <- "Lymphatic"
Idents(adata, cells = WhichCells(adata, idents = c("22"))) <- "Venule"
Idents(adata, cells = WhichCells(adata, idents = c("6", '17'))) <- "Capillary"
Idents(adata, cells = WhichCells(adata, idents = c("14"))) <- "Arteriole"

#annotate lymphocytes
Idents(adata, cells = WhichCells(adata, idents = c("33"))) <- "Proliferating T-cells"
Idents(adata, cells = WhichCells(adata, idents = c("2"))) <- "CD4"
Idents(adata, cells = WhichCells(adata, idents = c("3", "21"))) <- "CD8"
Idents(adata, cells = WhichCells(adata, idents = c("18", '27'))) <- "NKT"
Idents(adata, cells = WhichCells(adata, idents = c("7"))) <- "NK"
Idents(adata, cells = WhichCells(adata, idents = c("29"))) <- "B/Plasma"
Idents(adata, cells = WhichCells(adata, idents = c("31"))) <- "Mast"

#annotate myeloid cells
Idents(adata, cells = WhichCells(adata, idents = c("34"))) <- "Proliferating macrophage"
Idents(adata, cells = WhichCells(adata, idents = c("20"))) <- "moDC"
Idents(adata, cells = WhichCells(adata, idents = c("8", '10'))) <- "Monocyte-derived macrophage"
Idents(adata, cells = WhichCells(adata, idents = c("19", '15', '30', '11'))) <- "Monocyte"
Idents(adata, cells = WhichCells(adata, idents = c("16"))) <- "cDC"
Idents(adata, cells = WhichCells(adata, idents = c("5", '0', '1', '24_6'))) <- "Alveolar macrophage"

#save_annotated_file
adata$celltype <- Idents(adata)
saveRDS(adata, file = '~/Dropbox/Lab data/PCLS/pcls_merged_annotated_063021.rds')

#plot final annotation
DimPlot(adata, cols = 'polychrome', group.by = 'celltype')