# ---
# title: "CellLine_Atlas_Part1"
# author: "Zettabyte"
# ---


library(tidyverse)
library(ggplot2)
library(dplyr)
library(factoextra)
# library(RColorBrewer)
library(ggsci)
library(scales)
library(Rtsne)
library(stringr)
library(stats)
library(org.Hs.eg.db)
library(Seurat)
library(reshape2)
library(gtools)

set.seed(123)


season.colors1 <- c("#3c3241","#654a63","#6c5b74","#6d5f80","#8d7ba6","#9693af","#b6b3d5","#b1acbd","#c4c4b6","#eaead9")
barplot(rep(1,length(season.colors1)),col = season.colors1)
season.colors2 <- c("#1c1a19","#3b3636","#293425","#475a40","#a9c5b9","#c5dcd5","#dae9e7","#e6f0ef","#542447","#7a3467",
                    "#a15693","#c885c3","#ddb6e4","#fdcbff","#ffe0fe","#b06a6a","#d68181","#f3c2bd","#ecd6d9")
barplot(rep(1,length(season.colors2)),col = season.colors2)
season.colors3 <- c("#6c6f52","#91956e","#859374","#b6bb8a","#a7b992","#45767c","#5a9aa2","#6cb9c2","#37374d","#525273",
                    "#6d6d99","#9797d5","#41303a","#674c5c","#b694b4","#dcb3d9","#ffcff5","#b06a6a","#d68181","#f3c2bd","#ecd6d9")
barplot(rep(1,length(season.colors3)),col = season.colors3)

fig.colors <- c("#e64b35","#4dbbd6","#00a086","#3d5488","#f29b7f","#8491b4","#91d1c1","#dc0000","#7f6149","#F0C976",
                "#0044A2","#777F49","#b56f58","#00833e","#fbb045","#a5c2fc","#3a8ea3","#7f4019","#00d1bc","#ab65bc")
dark.fig.colors <- c("#b33a29","#3a8ea3","#006d5b","#263455","#bf7a64","#5e6781","#6d9e91","#a90000","#4c3a2b","#bd9e5c",
                     "#002e6e","#474c2b","#824f3f","#004f25","#c78c36","#839ac9","#27616f","#4c260e","#009e8e","#7c4989")
darker.fig.colors <- c("#80291d","#27616f","#003930","#0f1421","#8c5949","#383e4d","#496a62","#750000","#18120e","#8a7343",
                       "#00183b","#18180e","#4e3026","#001c0d","#936828","#617295","#15343c","#180b04","#006b60","#4d2d55")
# barplot(rep(1,length(fig.colors)),col = fig.colors)




setwd(output_dirname)



CL.choose <- "1"



#### 1.input ####
cellline <- paste0("xbx", CL.choose, "/")
cellline.project <- paste0("CL", CL.choose)

if (CL.choose=="3") {
  CL3.count <- read.csv("./XBX3WTA_RSEC_MolsPerCell.csv", header = TRUE,skip = 6)
  CL3.count <- CL3.count %>% 
    column_to_rownames("Cell_Index")
  CL <- CreateSeuratObject(counts = t(CL3.count), project = cellline.project, min.cell = 3, min.features = 200)
  CL
} else {
  CL <- CreateSeuratObject(Read10X(paste0(rawdata_dir, cellline)), project = cellline.project, min.cell = 3, min.features = 200)
  CL
}
if (CL.choose=="3") {
  CL[["percent.mt"]] <- PercentageFeatureSet(CL, pattern = "^MT\\.")
} else {
  CL[["percent.mt"]] <- PercentageFeatureSet(CL, pattern = "^MT-")
}
rb.genes <- rownames(CL)[grep("^RP[SL]",rownames(CL))]
CL.count <- GetAssayData(object = CL, slot = "counts")
percent.ribo <- Matrix::colSums(CL.count[rb.genes,])/Matrix::colSums(CL.count)*100
CL <- AddMetaData(CL, percent.ribo, col.name = "percent.ribo")

VlnPlot(CL, 
        features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), 
        group.by = "orig.ident", 
        cols = fig.colors, 
        pt.size = 0, ncol = 2) + labs(x="")
ggsave(paste0(cellline.project, "_1.QC_VlnPlot.pdf"), width = 5, height = 5)
ggsave(paste0(cellline.project, "_1.QC_VlnPlot.png"), width = 5, height = 5)


#### 2.Souporcell ####
spcell.CL <- read_tsv(paste0(souporcell_dir, cellline, "clusters_raw.tsv"))
spcell.cl2 <- read_tsv(paste0(souporcell_dir, cellline, "clusters_raw1000.tsv"))
spcell.cl3 <- read_tsv(paste0(souporcell_dir, cellline, "clusters_raw10000.tsv"))
spcell.cl4 <- read_tsv(paste0(souporcell_dir, cellline, "clusters_withKG.tsv"))
head(as.data.frame(spcell.CL), 2)

spcell.CL$cell_status <- paste0(spcell.CL$barcode,"_",spcell.CL$status)
spcell.cl2$cell_status <- paste0(spcell.cl2$barcode,"_",spcell.cl2$status)
spcell.cl3$cell_status <- paste0(spcell.cl3$barcode,"_",spcell.cl3$status)
spcell.cl4$cell_status <- paste0(spcell.cl4$barcode,"_",spcell.cl4$status)

spcell.CL$cell_assignment <- paste0(spcell.CL$barcode,"_",spcell.CL$assignment)
spcell.cl2$cell_assignment <- paste0(spcell.cl2$barcode,"_",spcell.cl2$assignment)
spcell.cl3$cell_assignment <- paste0(spcell.cl3$barcode,"_",spcell.cl3$assignment)
spcell.cl4$cell_assignment <- paste0(spcell.cl4$barcode,"_",spcell.cl4$assignment)
dim(spcell.CL)
test100 <- spcell.CL %>% filter(status=="singlet")
test100 <- split(test100$barcode, test100$assignment)
names(test100) <- paste0("iter100_", names(test100))
test1000 <- spcell.cl2 %>% filter(status=="singlet")
test1000 <- split(test1000$barcode, test1000$assignment)
test10000 <- spcell.cl3 %>% filter(status=="singlet")
test10000 <- split(test10000$barcode, test10000$assignment)
testKG <- spcell.cl4 %>% filter(status=="singlet")
testKG <- split(testKG$barcode, testKG$assignment)
names(testKG) <- paste0("KG_", names(testKG))

test <- lapply(lapply(testKG, function(x) lapply(test10000, function(y) {intersect(x, y) %>% length})), as.data.frame)
test.mtx <- tibble(test[[1]])
for (i in names(test)) {
  test.mtx <- rbind(test.mtx, test[[i]])
}
test.mtx <- test.mtx[-1,] %>% as.data.frame()
rownames(test.mtx) <- names(test)
test.mtx

test.mtx <- lapply(testKG, length) %>% as.data.frame() %>% t



#### 2.1.Cell status ####
library(VennDiagram)
venn.plot <- venn.diagram(
  x = list(Low_interation = spcell.CL$cell_status, Middle_interation = spcell.cl2$cell_status, 
           High_interation = spcell.cl3$cell_status, With_ref = spcell.cl4$cell_status),
  filename = paste0(output_dirname, cellline.project, "_2.VennDiagram.DiffPara_status.png"),width = 6000,height = 4500,
  lwd = 4,
  # col = "black",
  lty = 2,
  col = fig.colors[1:4],
  fill = fig.colors[1:4],
  alpha = 0.45,
  label.col = "black",
  cex = 2,
  fontfamily = "Arial",
  fontface = "bold",
  cat.col = dark.fig.colors[1:4],
  cat.cex = 2.2,
  cat.fontfamily = "Arial", 
  margin = 0.3 ,
  main = "Cell Line 1", 
  main.cex = 2, 
  main.pos = c(0.5, 0.95)
)


head(CL@meta.data)

CL$souporcell_status_100 <- NA
table(spcell.CL$status)
length(colnames(CL)[colnames(CL) %in% spcell.CL[spcell.CL$status%in%"doublet",]$barcode])
CL$souporcell_status_100[which(colnames(CL) %in% spcell.CL[spcell.CL$status%in%"doublet",]$barcode)] <- "Doublet"
length(colnames(CL)[colnames(CL) %in% spcell.CL[spcell.CL$status%in%"singlet",]$barcode])
CL$souporcell_status_100[which(colnames(CL) %in% spcell.CL[spcell.CL$status%in%"singlet",]$barcode)] <- "Singlet"
length(colnames(CL)[colnames(CL) %in% spcell.CL[spcell.CL$status%in%"unassigned",]$barcode])
CL$souporcell_status_100[which(colnames(CL) %in% spcell.CL[spcell.CL$status%in%"unassigned",]$barcode)] <- "Unassigned"
table(CL$souporcell_status_100)

CL$souporcell_status_1000 <- NA
table(spcell.cl2$status)
length(colnames(CL)[colnames(CL) %in% spcell.cl2[spcell.cl2$status%in%"doublet",]$barcode])
CL$souporcell_status_1000[which(colnames(CL) %in% spcell.cl2[spcell.cl2$status%in%"doublet",]$barcode)] <- "Doublet"
length(colnames(CL)[colnames(CL) %in% spcell.cl2[spcell.cl2$status%in%"singlet",]$barcode])
CL$souporcell_status_1000[which(colnames(CL) %in% spcell.cl2[spcell.cl2$status%in%"singlet",]$barcode)] <- "Singlet"
length(colnames(CL)[colnames(CL) %in% spcell.cl2[spcell.cl2$status%in%"unassigned",]$barcode])
CL$souporcell_status_1000[which(colnames(CL) %in% spcell.cl2[spcell.cl2$status%in%"unassigned",]$barcode)] <- "Unassigned"
table(CL$souporcell_status_1000)

CL$souporcell_status_10000 <- NA
table(spcell.cl3$status)
length(colnames(CL)[colnames(CL) %in% spcell.cl3[spcell.cl3$status%in%"doublet",]$barcode]) # 6508
CL$souporcell_status_10000[which(colnames(CL) %in% spcell.cl3[spcell.cl3$status%in%"doublet",]$barcode)] <- "Doublet"
length(colnames(CL)[colnames(CL) %in% spcell.cl3[spcell.cl3$status%in%"singlet",]$barcode]) # 27923
CL$souporcell_status_10000[which(colnames(CL) %in% spcell.cl3[spcell.cl3$status%in%"singlet",]$barcode)] <- "Singlet"
length(colnames(CL)[colnames(CL) %in% spcell.cl3[spcell.cl3$status%in%"unassigned",]$barcode]) # 393
CL$souporcell_status_10000[which(colnames(CL) %in% spcell.cl3[spcell.cl3$status%in%"unassigned",]$barcode)] <- "Unassigned"
table(CL$souporcell_status_10000)

CL$souporcell_status_withKG <- NA
table(spcell.cl4$status)
length(colnames(CL)[colnames(CL) %in% spcell.cl4[spcell.cl4$status%in%"doublet",]$barcode]) # 7623
CL$souporcell_status_withKG[which(colnames(CL) %in% spcell.cl4[spcell.cl4$status%in%"doublet",]$barcode)] <- "Doublet"
length(colnames(CL)[colnames(CL) %in% spcell.cl4[spcell.cl4$status%in%"singlet",]$barcode]) # 27000
CL$souporcell_status_withKG[which(colnames(CL) %in% spcell.cl4[spcell.cl4$status%in%"singlet",]$barcode)] <- "Singlet"
length(colnames(CL)[colnames(CL) %in% spcell.cl4[spcell.cl4$status%in%"unassigned",]$barcode]) # 201
CL$souporcell_status_withKG[which(colnames(CL) %in% spcell.cl4[spcell.cl4$status%in%"unassigned",]$barcode)] <- "Unassigned"
table(CL$souporcell_status_withKG)


################ assignment
head(CL@meta.data)
CL$souporcell_assign <- CL$souporcell_status_withKG

assign.num.ls <- spcell.cl4$assignment[!grepl("/", spcell.cl4$assignment)] %>% unique() %>% mixedsort()
assign.num.ls
for (i in assign.num.ls) {
  CL$souporcell_assign[which(colnames(CL) %in% spcell.cl4[spcell.cl4$assignment==i,]$barcode)] <- i
}

table(CL$souporcell_assign)
table(CL$souporcell_assign[which(colnames(CL) %in% spcell.cl4[spcell.cl4$status%in%"unassigned",]$barcode)])
CL$souporcell_assign[which(colnames(CL) %in% spcell.cl4[spcell.cl4$status%in%"unassigned",]$barcode)] <- "Unassigned"
table(CL$souporcell_assign)


test <- table(CL$souporcell_assign) %>% as.data.frame()
test$Var1 <- factor(test$Var1, levels = unique(CL$souporcell_assign) %>% mixedsort())
ggplot(test, aes(x=Var1, y=Freq, fill=Var1)) + 
  scale_fill_manual(values = c(fig.colors[1:(length(unique(CL$souporcell_assign))-2)],"grey","dark grey")) + 
  geom_bar(stat = 'identity') + 
  geom_text(aes(label = Freq), position = position_dodge(0.9), vjust = -0.8, size = 3.5) +
  theme_classic() + 
  theme(text = element_text(size = 15)) +
  ggtitle(paste0("CellLine", CL.choose)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(face = "plain", colour = "black")) +
  labs(x="", y="CellNum") +
  NoLegend()
ggsave(paste0(output_dirname, cellline.project, "_3.Barplot.souporcell_assign.CellNum.Raw.CL.pdf"), width = 9,height = 6)
ggsave(paste0(output_dirname, cellline.project, "_3.Barplot.souporcell_assign.CellNum.Raw.CL.png"), width = 9,height = 6)





#### 3.Singlets & subset #### 
CL$souporcell_status <- NA
CL@meta.data$souporcell_status[CL@meta.data$souporcell_status_100=="Singlet" & 
                                 CL@meta.data$souporcell_status_1000=="Singlet" & 
                                 CL@meta.data$souporcell_status_10000=="Singlet" & 
                                 CL@meta.data$souporcell_status_withKG=="Singlet"] <- "Singlet"
table(CL$souporcell_status)
Idents(CL) <- "souporcell_status"
CL_f <- subset(CL, ident="Singlet")
CL_f
table(CL_f$souporcell_assign)

test <- table(CL_f$souporcell_assign) %>% as.data.frame()
test$Var1 <- factor(test$Var1, levels = unique(CL_f$souporcell_assign) %>% mixedsort())
ggplot(test, aes(x=Var1, y=Freq, fill=Var1)) + 
  scale_fill_manual(values = fig.colors) + 
  geom_bar(stat = 'identity') + 
  geom_text(aes(label = Freq), position = position_dodge(0.9), vjust = -0.8, size = 3.5) +
  theme_classic() + 
  theme(text = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(face = "plain", colour = "black")) +
  ggtitle(paste0("CellLine", CL.choose)) + 
  labs(x="Singlets in Souporcell assignments", y="CellNum") +
  NoLegend()
ggsave(paste0(output_dirname, cellline.project, "_4.Barplot.souporcell_assign.CellNum.Singlets.CL_f.pdf"), width = 9,height = 6)
ggsave(paste0(output_dirname, cellline.project, "_4.Barplot.souporcell_assign.CellNum.Singlets.CL_f.png"), width = 9,height = 6)


CL_f <- subset(CL_f, subset = (nFeature_RNA >= 1500 & nFeature_RNA <= 8000 & percent.mt <= 30 & nCount_RNA >= 2000 & nCount_RNA < 30000))
CL_f
table(CL_f$souporcell_assign)


test <- table(CL_f$souporcell_assign) %>% as.data.frame()
test$Var1 <- factor(test$Var1, levels = unique(CL_f$souporcell_assign) %>% mixedsort())
ggplot(test, aes(x=Var1, y=Freq, fill=Var1)) + 
  scale_fill_manual(values = fig.colors) + 
  geom_bar(stat = 'identity') + 
  geom_text(aes(label = Freq), position = position_dodge(0.9), vjust = -0.8, size = 3.5) +
  theme_classic() + 
  theme(text = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(face = "plain", colour = "black")) +
  ggtitle(paste0("CellLine", CL.choose)) + 
  labs(x="After QC", y="CellNum") +
  NoLegend()
ggsave(paste0(output_dirname, cellline.project, "_5.Barplot.souporcell_assign.CellNum.QC.CL_f.pdf"), width = 9,height = 6)
ggsave(paste0(output_dirname, cellline.project, "_5.Barplot.souporcell_assign.CellNum.QC.CL_f.png"), width = 9,height = 6)


VlnPlot(CL_f, 
        features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), 
        group.by = "souporcell_assign", 
        cols = fig.colors, 
        pt.size = 0, ncol = 2) & labs(x="")
ggsave(paste0(cellline.project, "_6.QC_VlnPlot.pdf"), width = 7, height = 4.5)
ggsave(paste0(cellline.project, "_6.QC_VlnPlot.png"), width = 7, height = 4.5)




#### 4.Preprocessing & clustering ####

DefaultAssay(CL_f) <- "RNA"

CL_f <- NormalizeData(CL_f, normalization.method = "LogNormalize", scale.factor = 10000)
CL_f <- FindVariableFeatures(CL_f, selection.method = "vst", nfeatures = 5000)
CL_f <- ScaleData(CL_f, features = VariableFeatures(CL_f), verbose = TRUE)
CL_f <- RunPCA(CL_f, features = VariableFeatures(object = CL_f))
ifelse(CL.choose%in%c("1","2"), pick.res <- 0.2, pick.res <- 0.6)
ifelse(CL.choose=="1", pick.dims <- 10, pick.dims <- 20)

CL_f <- FindNeighbors(CL_f, dims = 1:pick.dims)
CL_f <- FindClusters(CL_f, resolution = pick.res)
CL_f <- RunUMAP(CL_f, dims = 1:pick.dims)
DimPlot(CL_f, label = T, repel = T, group.by = "seurat_clusters",reduction="umap",cols = c(fig.colors, "black","grey"))
ggsave(paste0(cellline.project, "_7.UMAP_CL_f_seurat_clusters.pdf"), width = 6, height = 4)
ggsave(paste0(cellline.project, "_7.UMAP_CL_f_seurat_clusters.png"), width = 6, height = 4)
DimPlot(CL_f, label = T, repel = T, group.by = "souporcell_assign",reduction="umap",cols = c(fig.colors, "black","grey"))
ggsave(paste0(cellline.project, "_7.UMAP_CL_f_souporcell_assign.pdf"), width = 6, height = 4)
ggsave(paste0(cellline.project, "_7.UMAP_CL_f_souporcell_assign.png"), width = 6, height = 4)

################ With Ref
head(CL_f@meta.data)
CL_f$souporcell_assign_WithRef <- CL_f$souporcell_status_withKG

assign.num.ls <- spcell.cl4$assignment[!grepl("/", spcell.cl4$assignment)] %>% unique() %>% mixedsort()
assign.num.ls
for (i in assign.num.ls) {
  CL_f$souporcell_assign_WithRef[which(colnames(CL_f) %in% spcell.cl4[spcell.cl4$assignment==i,]$barcode)] <- i
}

DimPlot(CL_f, label = T, repel = T, group.by = "souporcell_assign_WithRef",reduction="umap",cols = c(fig.colors, "black","grey"))
Idents(CL_f) <- "souporcell_assign_WithRef"
DimPlot(CL_f, label = F, repel = T, group.by = "souporcell_assign_WithRef",reduction="umap", 
        cells.highlight = WhichCells(object = CL_f, idents = "6")) + labs(title = "With Ref - Cluster6") + NoLegend()
DimPlot(CL_f, label = F, repel = T, group.by = "souporcell_assign_WithRef",reduction="umap", 
        cells.highlight = WhichCells(object = CL_f, idents = "8")) + labs(title = "With Ref - Cluster8") + NoLegend()
DimPlot(CL_f, label = F, repel = T, group.by = "souporcell_assign_WithRef",reduction="umap", 
        cells.highlight = WhichCells(object = CL_f, idents = "11")) + labs(title = "With Ref - Cluster11") + NoLegend()


################ 10000 iters
head(CL_f@meta.data)
CL_f$souporcell_assign_10000 <- CL_f$souporcell_status_10000

assign.num.ls <- spcell.cl3$assignment[!grepl("/", spcell.cl3$assignment)] %>% unique() %>% mixedsort()
assign.num.ls
for (i in assign.num.ls) {
  CL_f$souporcell_assign_10000[which(colnames(CL_f) %in% spcell.cl3[spcell.cl3$assignment==i,]$barcode)] <- i
}

Idents(CL_f) <- "souporcell_assign_10000"
DimPlot(CL_f, label = F, repel = T, group.by = "souporcell_assign_10000",reduction="umap", 
        cells.highlight = WhichCells(object = CL_f, idents = "0")) + labs(title = "No Ref (10000iters) - Cluster0") + NoLegend()
DimPlot(CL_f, label = F, repel = T, group.by = "souporcell_assign_10000",reduction="umap", 
        cells.highlight = WhichCells(object = CL_f, idents = "7")) + labs(title = "No Ref (10000iters) - Cluster7") + NoLegend()


#### 5.RNAref_DESeq2 ####

### top30
results.DESeq2 <- readRDS(paste0(DESeq2_dir, cellline, strsplit(cellline, "/")[[1]][1], "_results.DESeq2.rds"))
res.edgeR <- readRDS(paste0(edgeR_dir, cellline, strsplit(cellline, "/")[[1]][1], "_results.edgeR.rds"))

deg.lsls <- list()
for (j in 1:length(results.DESeq2)) {
  i=1
  deg.ls <- list()
  while (length(deg.ls) < 30) {
    deg.ls <- append(deg.ls, results.DESeq2[[j]][i])
    deg.ls <- append(deg.ls, res.edgeR[[j]][i])
    deg.ls <- unique(deg.ls)
    i = i + 1
  }
  deg.lsls[[names(results.DESeq2[j])]] <- unlist(deg.ls)
}

Idents(CL_f) <- "souporcell_assign"

for (i in 1:length(deg.lsls)) {
  CL.name <- names(deg.lsls[i])
  
  CL_f <- AddModuleScore(
    object = CL_f,
    features = list(deg.lsls[[i]]),
    ctrl = 10,
    name = paste0("Module.", CL.name, ".")
  )
  FeaturePlot(CL_f,features = paste0("Module.", CL.name, ".1"),label = T, order = T, min.cutoff = 0.1) + labs(x="", subtitle = paste0("CellLine", CL.choose))
  VlnPlot(object = CL_f, features = paste0("Module.", CL.name, ".1"), pt.size = 0, cols = fig.colors) + NoLegend() + labs(x="", subtitle = paste0("CellLine", CL.choose))
}


#### 6.AUCell ####

ifelse(CL.choose %in% c("1", "2"), Idents(CL_f) <- "RNA_snn_res.0.2", Idents(CL_f) <- "RNA_snn_res.0.6")

library(AUCell)
CL_f_rankings <- AUCell_buildRankings(CL_f@assays$RNA@data[rownames(CL_f),])

nrow(CL_f_rankings)
cells_AUC <- AUCell_calcAUC(geneSets = deg.lsls, 
                            rankings = CL_f_rankings, 
                            aucMaxRank = ceiling(nrow(CL_f_rankings)*0.05))
# dev.off()

for (i in names(cells_AUC)) {
  aucs <- as.numeric(getAUC(cells_AUC)[i, ])
  CL_f[[paste0("AUC_", i)]] <- aucs
  FeaturePlot(CL_f, features = paste0("AUC_", i),label = T, label.size = 6, order = T, label.color = "#c4c4b6") + labs(x="", subtitle = paste0("CellLine", CL.choose)) + viridis::scale_color_viridis(option="A")
}


#### 7.Annotation ####

CL1.name <- c("BT474","HCC1937","MCF10A","MCF7","MDAMB231","SKBR3","T47D","A549","MRC5","NCIH446","NCIH460","PC9")
CL2.name <- c("HELA","SIHA","786O","ACHN","HUH7","PLCPRF5","SKHEP1","EC109","KYSE150","KYSE510","TE1","ASPC1","CAPAN2","MIAPACA2","PANC1","SW1990","BCPAP","T24")
CL3.name <- c("C33A","AN3CA","HEC1A","HEC1B","HCT116","HT29","RKO","SW620","KYSE410","A2780","SKOV3","DU145","LNCaP","PC3","A431","CASKI",
              "DLD1","AGS","NCIN87","J82")

ifelse(CL.choose=="1", CL.name <- CL1.name, ifelse(CL.choose=="2", CL.name <- CL2.name, CL.name <- CL3.name))
zsc.mtx <- tibble(CL=CL.name)

for (i in 0:(length(CL.name)-1)) {
  data <- c()
  for (j in 1:length(CL.name)) {
    data <- append(data, CL_f[[paste0("Module.", CL.name[j], ".1")]][,1][CL_f@meta.data$seurat_clusters==i] %>% mean())
  }
  
  z_scores <- (data-mean(data))/sd(data)
  zsc.mtx[,paste0("cluster",i)] <- z_scores
  zsc.mtx[,paste0("raw.cluster",i)] <- data
}
zsc.mtx.raw <- zsc.mtx[,c(1,grep("^raw",colnames(zsc.mtx)))]
zsc.mtx <- zsc.mtx[,c(1,grep("^clus",colnames(zsc.mtx)))]
zsc.mtx
as.numeric(zsc.mtx[1,2:ncol(zsc.mtx)]) %>% mean
as.numeric(zsc.mtx[1,2:ncol(zsc.mtx)]) %>% sd
zsc.mtx$cluster0 %>% mean
zsc.mtx$cluster0 %>% sd

test <- melt(zsc.mtx, id.vars = c("CL"))
test$CL <- factor(test$CL, levels = CL.name)
test <- arrange(test, CL)


library(igraph)
test1 <- test
test1 <- test1[test1$value > 0,]
test1 <- arrange(test1, desc(value))

p <- graph_from_data_frame(test1, directed = F)
# E(p)$width <- test1$value*3
V(p)$type <- as.logical(c(rep(0, length(CL.name)), rep(1, length(CL.name))))

netm <- get.adjacency(p, attr="value", sparse=F)
netm <- netm[1:length(CL.name), (length(CL.name)+1):ncol(netm)]
netm <- netm[,mixedsort(colnames(netm))]
netm <- netm[CL.name,]
netm <- netm[order(apply(netm, 1, max)), order(apply(netm, 2, max))]
write.csv(netm, paste0(output_dirname, cellline.project, "_10.iGraph_adjacency.csv"))

ggplot(melt(netm), aes(x=Var2, y=Var1, fill=value)) +
  geom_tile() +
  labs(x="Column", y="Row") +
  scale_fill_gradient(low=fig.colors[2], high=fig.colors[1]) + # 调整颜色渐变范围
  theme_void() +
  labs(x="", y="") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = element_text(size = 12, family = "", face = "plain", colour = "black"),
        plot.margin = unit(rep(1,4),"cm"),
        legend.position = "none")


##### VlnPlot MGS & AUC featureplot #####

test.idx <- grep("Module", colnames(CL1_f@meta.data))
colnames(CL1_f@meta.data)[test.idx]
test <- gsub(pattern = "Module.", replacement = "MGS_", x = colnames(CL1_f@meta.data)[test.idx])
colnames(CL1_f@meta.data)[test.idx] <- gsub(pattern = "\\.1", replacement = "", x = test)

CL1_f$RNA_snn_res.0.2 <- factor(CL1_f$RNA_snn_res.0.2, levels = c(6, 8, 10, 7, 1, 11, 9, 2, 4, 3, 5, 0))
VlnPlot(CL1_f,features = grep("MGS", colnames(CL1_f@meta.data), value = T), 
        stack = T, flip = T, group.by = "RNA_snn_res.0.2", cols = c(fig.colors, dark.fig.colors)) + 
  labs(x="") + NoLegend()


test.idx <- grep("Module", colnames(CL2_f@meta.data))
colnames(CL2_f@meta.data)[test.idx]
test <- gsub(pattern = "Module.", replacement = "MGS_", x = colnames(CL2_f@meta.data)[test.idx])
colnames(CL2_f@meta.data)[test.idx] <- gsub(pattern = "\\.1", replacement = "", x = test)

CL2_f$RNA_snn_res.0.2 <- factor(CL2_f$RNA_snn_res.0.2, levels = c(11, 13, 17, 7, 4, 1, 2, 14, 10, 12, 5, 0, 3, 16, 9, 15, 6, 8))
VlnPlot(CL2_f,features = grep("MGS", colnames(CL2_f@meta.data), value = T), 
        stack = T, flip = T, group.by = "RNA_snn_res.0.2", cols = c(fig.colors, dark.fig.colors)) + 
  labs(x="") + NoLegend()



test.idx <- grep("Module", colnames(CL3_f@meta.data))
colnames(CL3_f@meta.data)[test.idx]
test <- gsub(pattern = "Module.", replacement = "MGS_", x = colnames(CL3_f@meta.data)[test.idx])
colnames(CL3_f@meta.data)[test.idx] <- gsub(pattern = ".1", replacement = "", x = test)

CL3_f$RNA_snn_res.0.6 <- factor(CL3_f$RNA_snn_res.0.6, levels = c(9, 1, 10, 8, 16, 17, 7, 0, 4, 3, 19, 6, 2, 15, 11, 14, 13, 18, 5, 12))
VlnPlot(CL3_f,features = grep("MGS", colnames(CL3_f@meta.data), value = T), 
        stack = T, flip = T, group.by = "RNA_snn_res.0.6", cols = c(fig.colors, dark.fig.colors)) + 
  labs(x="") + NoLegend()




