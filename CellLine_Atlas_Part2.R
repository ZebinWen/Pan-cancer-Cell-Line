# ---
# title: "CellLine_Atlas_Part2"
# author: "Zettabyte"
# cite: "Pan-cancer single-cell RNA-seq identifies recurring programs of cellular heterogeneity" and "Hallmarks of transcriptional intratumour heterogeneity across a thousand tumours"
# ---




library(tidyverse)
library(ggplot2)
library(dplyr)
library(factoextra)
library(ggsci)
library(scales)
library(Rtsne)
library(stringr)
library(stats)
library(org.Hs.eg.db)
library(Seurat)
library(reshape2)
library(gtools)
library(RColorBrewer)
library(viridis)
library(enrichplot)
library(clusterProfiler)
library(aplot)

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

library(paletteer)
my_color = paletteer_d(`"ggsci::default_nejm"`)
my_color = colorRampPalette(my_color)(50)
my_color


setwd(output_dirname)


theme1 <- function(..., bg="white"){
  require(grid)
  theme_classic(...) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = bg, colour = "black", size = 1.2),
      legend.position = "right",
      legend.text = element_text(size = 12),
      legend.key = element_rect(fill = bg, color = NA),
      legend.title = element_text(color = "black", size = 13),
      axis.text.x = element_text(angle = 0, colour = "black", size = 12),
      axis.text.y = element_text(angle = 0, colour = "black", size = 12),
      plot.title = element_text(hjust = 0.5, size = 14, colour = "black"),
      axis.title.x = element_text(size = 12, colour = "black"),
      axis.title.y = element_text(size = 12, colour = "black")
    )
}




#### 1.adjacency/SNPvsRNA heatmap ####

netm <- read.csv("./CL_ALL_10.iGraph_adjacency.csv")
netm <- column_to_rownames(netm, "X")
head(netm)
apply(netm, 2, max) %>% sort(decreasing = T)
netm <- netm[names(apply(netm, 1, max) %>% sort(decreasing = F)),]
netm <- netm[,names(apply(netm, 2, max) %>% sort(decreasing = F))]
netm$X <- rownames(netm)
test <- melt(netm)
test$X <- factor(test$X, levels = netm$X)

test <- test[!(test$X %in% c("EC109", "PLCPRF5", "HELA", "KYSE510", "KYSE410")),]
test <- test[!(test$variable %in% c("CL2_cluster0", "CL2_cluster1", "CL2_cluster3", "CL2_cluster10", "CL3_cluster11")),]

ggplot(test, aes(x=variable, y=X, fill=value)) +
  geom_tile() +
  labs(x="Column", y="Row") +
  scale_fill_gradient(low=fig.colors[2], high=fig.colors[1]) +
  theme_void() +
  labs(x="", y="") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = element_text(size = 12, family = "", face = "plain", colour = "black"),
        plot.margin = unit(rep(1,4),"cm"),
        legend.position = "none")

CL1.df <- test[str_detect(test$variable, "CL1"),]
CL1.df <- CL1.df %>% group_by(variable) %>% top_n(n=1, wt=value)
CL1.df$seurat_cluster <- lapply(CL1.df$variable, function(x){strsplit(x %>% as.character(), split = "er")[[1]][2]}) %>% unlist()
CL1.df
CL2.df <- test[str_detect(test$variable, "CL2"),]
CL2.df <- CL2.df %>% group_by(variable) %>% top_n(n=1, wt=value)
CL2.df$seurat_cluster <- lapply(CL2.df$variable, function(x){strsplit(x %>% as.character(), split = "er")[[1]][2]}) %>% unlist()
CL2.df
CL3.df <- test[str_detect(test$variable, "CL3"),]
CL3.df <- CL3.df %>% group_by(variable) %>% top_n(n=1, wt=value)
CL3.df$seurat_cluster <- lapply(CL3.df$variable, function(x){strsplit(x %>% as.character(), split = "er")[[1]][2]}) %>% unlist()
CL3.df



#### matching heatmap
test <- read.csv("6.SNP_RNA_assign_results_CL1.csv")
test1 <- column_to_rownames(test, "X")
test1 <- test1 / rowSums(test1)
test <- cbind(test[,1], test1)
colnames(test)[1] <- "RNA_CellLine"
test.cl1 <- melt(test)
test.cl1$variable <- paste0("CL1_SNP_", test.cl1$variable)
test.cl1$RNA_CellLine <- paste0("CL1_", test.cl1$RNA_CellLine)

test <- read.csv("6.SNP_RNA_assign_results_CL2.csv")
test1 <- column_to_rownames(test, "X")
test1 <- test1 / rowSums(test1)
test <- cbind(test[,1], test1)
colnames(test)[1] <- "RNA_CellLine"
colnames(test)[2] <- "786O"
colnames(test)[2:ncol(test)] <- paste0("CL2_SNP_", colnames(test)[2:ncol(test)])
test.cl2 <- melt(test)
test.cl2$RNA_CellLine <- paste0("CL2_", test.cl2$RNA_CellLine)

test <- read.csv("6.SNP_RNA_assign_results_CL3.csv")
test1 <- column_to_rownames(test, "X")
test1 <- test1 / rowSums(test1)
test <- cbind(test[,1], test1)
colnames(test)[1] <- "RNA_CellLine"
test.cl3 <- melt(test)
test.cl3$variable <- paste0("CL3_SNP_", test.cl3$variable)
test.cl3$RNA_CellLine <- paste0("CL3_", test.cl3$RNA_CellLine)

test <- rbind(test.cl1, test.cl2)
test <- rbind(test, test.cl3)

test$RNA_CellLine <- lapply(test$RNA_CellLine, function(x){strsplit(x, "_")[[1]][3]}) %>% unlist
test$RNA_CellLine <- factor(test$RNA_CellLine, levels = test$RNA_CellLine[c(1:12, 145:158, 341:357)])
test$variable <- lapply(test$variable, function(x){strsplit(x, "_")[[1]][3]}) %>% unlist
test$variable <- factor(test$variable, levels = test$RNA_CellLine[c(1:12, 145:158, 341:357)])
test

ggplot(test, aes(x=variable, y=RNA_CellLine, fill=value)) +
  geom_tile() +
  geom_rect(aes(xmin = 0.5, xmax = 12.5, ymin = 0.5, ymax = 12.5), color = "black", fill = NA) +
  geom_rect(aes(xmin = 12.5, xmax = 26.5, ymin = 12.5, ymax = 26.5), color = "black", fill = NA) +
  geom_rect(aes(xmin = 26.5, xmax = 43.5, ymin = 26.5, ymax = 43.5), color = "black", fill = NA) +
  geom_text(aes( x=3, y=16, label="CL1"), color="black", size=4.5 , angle=0, fontface="plain") +
  geom_text(aes( x=14, y=29.5, label="CL2"), color="black", size=4.5 , angle=0, fontface="plain") +
  geom_text(aes( x=24, y=40, label="CL3"), color="black", size=4.5 , angle=0, fontface="plain") +
  scale_fill_gradient(low=fig.colors[2], high=fig.colors[1]) +
  theme_bw() +
  labs(x="", y="", title = "Demultiplexed Results (SNP vs RNA)", fill="") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = element_text(size = 12, family = "", face = "plain", colour = "black"),
        plot.margin = unit(rep(1,4),"cm"),
        )

ggplot(test, aes(x=variable, y=RNA_CellLine, fill=value)) +
  geom_tile() +
  geom_rect(aes(xmin = 0.5, xmax = 12.5, ymin = 0.5, ymax = 12.5), color = "black", fill = NA) +
  geom_rect(aes(xmin = 12.5, xmax = 26.5, ymin = 12.5, ymax = 26.5), color = "black", fill = NA) +
  geom_rect(aes(xmin = 26.5, xmax = 43.5, ymin = 26.5, ymax = 43.5), color = "black", fill = NA) +
  geom_text(aes( x=3, y=16, label="CL1"), color="black", size=4.5 , angle=0, fontface="plain") +
  geom_text(aes( x=14, y=29.5, label="CL2"), color="black", size=4.5 , angle=0, fontface="plain") +
  geom_text(aes( x=24, y=40, label="CL3"), color="black", size=4.5 , angle=0, fontface="plain") +
  scale_fill_gradient(low=fig.colors[2], high=fig.colors[1]) + 
  theme_bw() +
  labs(x="", y="", title = "Demultiplexed Results (SNP vs RNA)", fill="") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
        axis.ticks.length.y = unit(x = 0, units = "cm"),
        axis.text = element_text(size = 12, family = "", face = "plain", colour = "black"),
        plot.margin = unit(rep(1,4),"cm"),
  )


#### 2.Confirm Celltype ####

CL1.name <- c("BT474","HCC1937","MCF10A","MCF7","MDAMB231","SKBR3","T47D","A549","MRC5","NCIH446","NCIH460","PC9")
CL2.name <- c("HELA","SIHA","786O","ACHN","HUH7","PLCPRF5","SKHEP1","EC109","KYSE150","KYSE510","TE1","ASPC1","CAPAN2","MIAPACA2","PANC1","SW1990","BCPAP","T24")
CL3.name <- c("C33A","AN3CA","HEC1A","HEC1B","HCT116","HT29","RKO","SW620","KYSE410","A2780","SKOV3","DU145","LNCaP","PC3","A431","CASKI",
              "DLD1","AGS","NCIN87","J82")

CL1_f <- readRDS("./CL1_f.rds")
DimPlot(CL1_f, group.by = "souporcell_CellLine", cols = fig.colors) + NoLegend() | DimPlot(CL1_f, group.by = "RNA_CellLine", cols = fig.colors)
CL2_f <- readRDS("./CL2_f.rds")
DimPlot(CL2_f, group.by = "souporcell_CellLine", cols = fig.colors) + NoLegend() | DimPlot(CL2_f, group.by = "RNA_CellLine", cols = fig.colors)
CL3_f <- readRDS("./CL3_f.rds")
DimPlot(CL3_f, group.by = "souporcell_CellLine", cols = fig.colors) + NoLegend() | 
  DimPlot(CL3_f, group.by = "RNA_CellLine", cols = fig.colors)

CL1_f$CellLine <- CL1_f$RNA_CellLine
CL2_f$CellLine <- CL2_f$RNA_CellLine
CL3_f$CellLine <- CL3_f$RNA_CellLine
CL3_f$CellLine[CL3_f$CellLine=="LNCaP"] <- "LNCAP"

CL1_f <- RenameCells(CL1_f, add.cell.id = "CL1")
CL1_f$disease <- "Cancer"
CL1_f$disease[CL1_f$CellLine%in%c("MRC5", "MCF10A")] <- "Normal"
CL2_f <- RenameCells(CL2_f, add.cell.id = "CL2")
CL2_f$disease <- "Cancer"
CL3_f <- RenameCells(CL3_f, add.cell.id = "CL3")
CL3_f$disease <- "Cancer"



#### 3.Merge to one object ####

CL <- merge(CL1_f, c(CL2_f, CL3_f))

DefaultAssay(CL) <- "RNA"
CL <- NormalizeData(CL, normalization.method = "LogNormalize", scale.factor = 10000)
CL <- FindVariableFeatures(CL, selection.method = "vst", nfeatures = 7000)
CL <- ScaleData(CL, features = VariableFeatures(CL), verbose = TRUE)
CL <- RunPCA(CL, features = VariableFeatures(object = CL))
CL <- FindNeighbors(CL, dims = 1:20)
CL <- FindClusters(CL, resolution = 0.8)
CL <- RunUMAP(CL, dims = 1:20)
CL <- RunTSNE(CL, dims = 1:20)
DimPlot(CL, label = T, repel = T, group.by = "CellLine",reduction="umap",cols = c(fig.colors, dark.fig.colors, darker.fig.colors)) + NoLegend()

Idents(CL) <- "CellLine"
umap_coords <- as.data.frame(CL@reductions$umap@cell.embeddings)
colnames(umap_coords) <- c("UMAP1", "UMAP2")
umap_coords$cluster <- as.character(Idents(CL))
test <- umap_coords %>% group_by(cluster) %>% summarize(avg_1=mean(UMAP1), avg_2=mean(UMAP2))

library(ggrepel)
ggplot(umap_coords, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = c(fig.colors, dark.fig.colors, darker.fig.colors)) +
  theme_bw() +
  geom_text_repel(data = test, aes(x=avg_1, y=avg_2, label=cluster), color = "black", bg.color = "white", bg.r = 0.15) +
  theme(legend.position = "None")




#### 4.Cell line info ####

CL$tissue <- NA
CL$tissue[CL$CellLine%in%c("A549", "PC9", "NCIH446", "NCIH460", "MRC5")] <- "Lung Cancer"
CL$tissue[CL$CellLine%in%c("SKBR3", "HCC1937", "MCF10A", "BT474", "T47D", "MDAMB231", "MCF7")] <- "Breast Cancer"
CL$tissue[CL$CellLine%in%c("HUH7", "SKHEP1", "PLCPRF5")] <- "Liver Cancer"
CL$tissue[CL$CellLine%in%c("TE1", "EC109", "KYSE510", "KYSE150",
                           "KYSE410")] <- "Esophageal Cancer"
CL$tissue[CL$CellLine%in%c("HELA", "SIHA",
                           "C33A", "CASKI")] <- "Cervical Cancer"
CL$tissue[CL$CellLine%in%c("BCPAP")] <- "Thyroid Cancer"
CL$tissue[CL$CellLine%in%c("T24", 
                           "J82")] <- "Bladder Cancer"
CL$tissue[CL$CellLine%in%c("786O", "ACHN")] <- "Kidney Cancer"
CL$tissue[CL$CellLine%in%c("PANC1", "ASPC1", "SW1990", "CAPAN2", "MIAPACA2")] <- "Pancreatic Cancer"
CL$tissue[CL$CellLine%in%c("DU145", "PC3", "LNCAP")] <- "Prostate Cancer"
CL$tissue[CL$CellLine%in%c("A2780", "SKOV3")] <- "Ovarian Cancer"
CL$tissue[CL$CellLine%in%c("AGS", "NCIN87")] <- "Gastric Cancer"
CL$tissue[CL$CellLine%in%c("A431")] <- "Skin Cancer"
CL$tissue[CL$CellLine%in%c("HEC1B", "HEC1A", "AN3CA")] <- "Endometrial/Uterine Cancer"
CL$tissue[CL$CellLine%in%c("RKO", "SW620", "HT29", "HCT116", "DLD1")] <- "Colon/Colorectal Cancer"

CL$disease <- "Cancer"
CL$disease[CL$CellLine%in%c("MRC5", "MCF10A")] <- "Normal"

DimPlot(CL, group.by = "tissue", label = T, repel = T, cols = c(season.colors1, season.colors2))
DimPlot(CL, group.by = "tissue", label = F, cols = fig.colors) + labs(title = "Cancer Type")
DimPlot(CL, group.by = "CellLine", label = T, repel = T, cols = c(fig.colors, dark.fig.colors, darker.fig.colors))




#### 5.NMF ####
## (In this section, some codes were cited from the paper "Pan-cancer single-cell RNA-seq identifies recurring programs of cellular heterogeneity" and "Hallmarks of transcriptional intratumour heterogeneity across a thousand tumours")
rm_lowCL <- function(obj, min_num=100, obj.ident="CellLine") {
  target.CLs <- table(obj[[obj.ident]])[table(obj[[obj.ident]]) > 100] %>% names()
  Idents(obj) <- obj.ident
  obj2 <- subset(obj, idents=target.CLs)
  return(obj2)
}
get_expr_list <- function(obj, obj.ident="CellLine", output.name=NULL) {
  if (is.null(output.name)) {
    stop("Error: Please type the output name!")
  }
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 7000)
  target.genes <- VariableFeatures(obj)
  Idents(obj) <- obj.ident
  expr_list <- list()
  for (i in table(obj[[obj.ident]]) %>% names()) {
    print(paste0("#### ", i , " ####"))
    test <- subset(obj, ident=i)
    test.count <- as.data.frame(test@assays$RNA@counts)
    mtx <- apply(test.count[target.genes, ], 2, function(x) {
      x / sum(x) * 1000000
    })
    mtx <- log2((mtx/10)+1)
    print(paste0("Gene numbers before removed: ", dim(mtx)[1]))
    mtx <- mtx[apply(mtx, 1, function(x) { length(which(x > 3.5)) > ncol(mtx) * 0.02 }), ]
    print(paste0("Gene numbers after removed: ", dim(mtx)[1]))
    mtx <- mtx - rowMeans(mtx)
    mtx[mtx < 0] = 0
    expr_list[[i]] <- mtx
  }
  saveRDS(expr_list, output.name)
}




#### 5.1 CL原始矩阵预处理 ####

for (cl in c("CL1_f", "CL2_f", "CL3_f")) {
  CL <- get(cl)
  CL <- FindVariableFeatures(CL, selection.method = "vst", nfeatures = 7000)
  target.genes <- VariableFeatures(CL)
  
  Idents(CL) <- "CellLine"
  expr_list <- list()
  for (i in (table(CL$CellLine) %>% names())) {
    # i = "A2780"
    print(paste0("#### ", i , " ####"))
    test <- subset(CL, ident=i)
    test.count <- as.data.frame(test@assays$RNA@counts)
    mtx <- apply(test.count[target.genes, ], 2, function(x) {
      x / sum(x) * 1000000
    })
    mtx <- log2((mtx/10)+1)
    mtx <- mtx[apply(mtx, 1, function(x) { length(which(x > 3.5)) > ncol(mtx) * 0.02 }), ]
    print(paste0("Gene numbers after removed: ", dim(mtx)[1]))
    mtx <- mtx - rowMeans(mtx)
    mtx[mtx < 0] = 0
    expr_list[[i]] <- mtx
  }
  saveRDS(expr_list, paste0("4.", cl, "_expr_list.rds"))
}



#### 5.2 Public data ####

library(Matrix)
library(data.table)
library(AnnoProbe)


#### 5.2.1 BC32 - 2022 ####

BC32 <- readRDS(paste0(public_dir, "GSE173634/RAW.UMI.counts.BC.cell.lines.rds"))
ids=annoGene(rownames(BC32),'ENSEMBL','human')
head(ids)
ids=ids[!duplicated(ids$SYMBOL),]
ids=ids[!duplicated(ids$ENSEMBL),]
dim(ids)
kp = rownames(BC32) %in% ids$ENSEMBL
table(kp)
BC32 <- BC32[kp,]
rownames(BC32)=ids$SYMBOL[match( rownames(BC32) , ids$ENSEMBL)]
BC32[1:4,1:4]
BC32[1:4,(ncol(BC32)-4):ncol(BC32)]

BC32.seu <- CreateSeuratObject(counts = BC32, min.cells = 3, min.features = 200)
BC32.seu$CellLine <- lapply(colnames(BC32.seu), function(x){
  strsplit(x, "_")[[1]][1]
}) %>% unlist()


table(BC32.seu[["CellLine"]])[table(BC32.seu[["CellLine"]]) > 100] %>% names() %>% length()

get_expr_list(obj = BC32.seu, obj.ident = "CellLine", output.name = paste0("4.", "BC32.seu", "_expr_list.rds"))





#### 5.2.2 BC6 - 2021 ####

BC6 <- read.csv(paste0(public_dir, "GSM4285803/GSM4285803_scRNA_RawCounts.csv"))
BC6 <- column_to_rownames(BC6, "X")
BC6 <- t(BC6)
BC6[1:4,1:4]
colnames(BC6)
dim(BC6) # 21888 4614

mtdt <- read.csv(paste0(public_dir, "GSM4285803/GSM4285803_scRNA_metaInfo.csv"))
mtdt <- column_to_rownames(mtdt, "X")
head(mtdt);dim(mtdt)

BC6.seu <- CreateSeuratObject(BC6, project = "BC6_2021", meta.data = mtdt, min.cells = 3, min.features = 200)
BC6.seu$CellLine <- BC6.seu$CellType
table(BC6.seu$CellLine)
Idents(BC6.seu) <- "CellLine"
BC7.seu <- subset(BC6.seu, ident="T47D KO", invert=T)
BC7.seu$CellLine[BC7.seu$CellLine=="T47D WT"] <- "T47D"

get_expr_list(obj = BC7.seu, obj.ident = "CellLine", output.name = paste0("4.", "BC7.seu", "_expr_list.rds"))






#### 5.2.3 CCLE - 2020 ####
CCLE <- read.table(paste0(public_dir, "CCLE/UMIcount_data.txt"))
CCLE[1:4,1:4]

mtdt <- CCLE[1:2,]
mtdt[1:2,1:4]
mtdt <- t(mtdt) %>% as.data.frame()
unique(mtdt$Cell_line); dim(mtdt) # 56982     2

CCLE <- CCLE[-c(1:2),]
dim(CCLE) # 30314 56982

CCLE.seu <- CreateSeuratObject(CCLE, project = "CCLE_2020", meta.data = mtdt, min.cells = 3, min.features = 200)
CCLE.seu$CellLine <- CCLE.seu$Cell_line
CCLE.seu@meta.data$Cell_line <- lapply(CCLE.seu@meta.data$Cell_line, function(x){strsplit(x, split = "_")[[1]][1]}) %>% unlist()

CCLE.seu[["percent.mt"]] <- PercentageFeatureSet(CCLE.seu, pattern = "^MT-")
CCLE.seu <- subset(CCLE.seu, subset = (nFeature_RNA >= 1500 & nFeature_RNA <= 8000 & percent.mt <= 30 & nCount_RNA >= 2000 & nCount_RNA < 30000))
CCLE.seu <- NormalizeData(CCLE.seu, normalization.method = "LogNormalize", scale.factor = 10000)
CCLE.seu <- FindVariableFeatures(CCLE.seu, selection.method = "vst", nfeatures = 5000)
CCLE.seu <- ScaleData(CCLE.seu, features = VariableFeatures(CCLE.seu), verbose = TRUE)
CCLE.seu <- RunPCA(CCLE.seu, features = VariableFeatures(object = CCLE.seu))
CCLE.seu <- FindNeighbors(CCLE.seu, dims = 1:20)
CCLE.seu <- FindClusters(CCLE.seu, resolution = 0.8)
CCLE.seu <- RunUMAP(CCLE.seu, dims = 1:20)
DimPlot(CCLE.seu, label = T, repel = T, group.by = "seurat_clusters",reduction="umap",cols = c(fig.colors, "black","grey"))
ggsave(paste0("UMAP_CCLE.seu_seurat_clusters.pdf"), width = 6, height = 4)
ggsave(paste0("UMAP_CCLE.seu_seurat_clusters.png"), width = 6, height = 4)


Idents(CCLE.seu) <- "Pool_ID"
CCLE.seu.Pool18 <- subset(CCLE.seu, ident="18")
CCLE.seu.Pool16 <- subset(CCLE.seu, ident="16")
CCLE.seu.Poolcustom <- subset(CCLE.seu, ident="custom")
CCLE.seu.Pool22 <- subset(CCLE.seu, ident="22")
CCLE.seu.Pool19 <- subset(CCLE.seu, ident="19")
CCLE.seu.Pool9 <- subset(CCLE.seu, ident="9")
CCLE.seu.Pool15 <- subset(CCLE.seu, ident="15")
CCLE.seu.Pool10 <- subset(CCLE.seu, ident="10")
CCLE.seu.Pool6 <- subset(CCLE.seu, ident="6")

table(CCLE.seu.Pool18[["CellLine"]])[table(CCLE.seu.Pool18[["CellLine"]]) > 100] %>% names() %>% length()
CCLE.seu.Pool18.20 <- rm_lowCL(obj = CCLE.seu.Pool18, min_num = 100, obj.ident = "CellLine")

table(CCLE.seu.Pool16[["CellLine"]])[table(CCLE.seu.Pool16[["CellLine"]]) > 100] %>% names() %>% length()
CCLE.seu.Pool16.16 <- rm_lowCL(obj = CCLE.seu.Pool16, min_num = 100, obj.ident = "CellLine")

table(CCLE.seu.Poolcustom[["CellLine"]])[table(CCLE.seu.Poolcustom[["CellLine"]]) > 100] %>% names() %>% length()
CCLE.seu.Poolcustom.8 <- rm_lowCL(obj = CCLE.seu.Poolcustom, min_num = 100, obj.ident = "CellLine")

table(CCLE.seu.Pool22[["CellLine"]])[table(CCLE.seu.Pool22[["CellLine"]]) > 100] %>% names() %>% length()
CCLE.seu.Pool22.23 <- rm_lowCL(obj = CCLE.seu.Pool22, min_num = 100, obj.ident = "CellLine")

table(CCLE.seu.Pool19[["CellLine"]])[table(CCLE.seu.Pool19[["CellLine"]]) > 100] %>% names() %>% length()
CCLE.seu.Pool19.23 <- rm_lowCL(obj = CCLE.seu.Pool19, min_num = 100, obj.ident = "CellLine")

table(CCLE.seu.Pool9[["CellLine"]])[table(CCLE.seu.Pool9[["CellLine"]]) > 100] %>% names() %>% length()
CCLE.seu.Pool9.25 <- rm_lowCL(obj = CCLE.seu.Pool9, min_num = 100, obj.ident = "CellLine")

table(CCLE.seu.Pool15[["CellLine"]])[table(CCLE.seu.Pool15[["CellLine"]]) > 100] %>% names() %>% length()
CCLE.seu.Pool15.23 <- rm_lowCL(obj = CCLE.seu.Pool15, min_num = 100, obj.ident = "CellLine")

table(CCLE.seu.Pool10[["CellLine"]])[table(CCLE.seu.Pool10[["CellLine"]]) > 100] %>% names() %>% length()
CCLE.seu.Pool10.25 <- rm_lowCL(obj = CCLE.seu.Pool10, min_num = 100, obj.ident = "CellLine")

table(CCLE.seu.Pool6[["CellLine"]])[table(CCLE.seu.Pool6[["CellLine"]]) > 100] %>% names() %>% length()
CCLE.seu.Pool6.24 <- rm_lowCL(obj = CCLE.seu.Pool6, min_num = 100, obj.ident = "CellLine")


CCLE.seu.Pool18.20
get_expr_list(obj = CCLE.seu.Pool18.20, obj.ident = "CellLine", output.name = paste0("4.", "CCLE.seu.Pool18.20", "_expr_list.rds"))
CCLE.seu.Pool16.16
get_expr_list(obj = CCLE.seu.Pool16.16, obj.ident = "CellLine", output.name = paste0("4.", "CCLE.seu.Pool16.16", "_expr_list.rds"))
CCLE.seu.Poolcustom.8
get_expr_list(obj = CCLE.seu.Poolcustom.8, obj.ident = "CellLine", output.name = paste0("4.", "CCLE.seu.Poolcustom.8", "_expr_list.rds"))
CCLE.seu.Pool22.23
get_expr_list(obj = CCLE.seu.Pool22.23, obj.ident = "CellLine", output.name = paste0("4.", "CCLE.seu.Pool22.23", "_expr_list.rds"))
CCLE.seu.Pool19.23
get_expr_list(obj = CCLE.seu.Pool19.23, obj.ident = "CellLine", output.name = paste0("4.", "CCLE.seu.Pool19.23", "_expr_list.rds"))
CCLE.seu.Pool9.25
get_expr_list(obj = CCLE.seu.Pool9.25, obj.ident = "CellLine", output.name = paste0("4.", "CCLE.seu.Pool9.25", "_expr_list.rds"))
CCLE.seu.Pool15.23
get_expr_list(obj = CCLE.seu.Pool15.23, obj.ident = "CellLine", output.name = paste0("4.", "CCLE.seu.Pool15.23", "_expr_list.rds"))
CCLE.seu.Pool10.25
get_expr_list(obj = CCLE.seu.Pool10.25, obj.ident = "CellLine", output.name = paste0("4.", "CCLE.seu.Pool10.25", "_expr_list.rds"))
CCLE.seu.Pool6.24
get_expr_list(obj = CCLE.seu.Pool6.24, obj.ident = "CellLine", output.name = paste0("4.", "CCLE.seu.Pool6.24", "_expr_list.rds"))




#### 5.2.4 LCL - 2021 ####
LCL_777_B958 <- Read10X(data.dir = "./GSE158275_RAW_LCL/LCL_777_B958/")
LCL_777_B958 <- CreateSeuratObject(counts = LCL_777_B958, project = "LCL_777_B958", min.cells = 3, min.features = 200)
LCL_777_B958[["percent.mt"]] <- PercentageFeatureSet(LCL_777_B958, pattern = "^MT-")
VlnPlot(LCL_777_B958, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(LCL_777_B958, features = c("nFeature_RNA"), y.max = 2000)
VlnPlot(LCL_777_B958, features = c("nCount_RNA"), y.max = 8000)
VlnPlot(LCL_777_B958, features = c("percent.mt"), y.max = 10)

LCL_777_B958 <- subset(LCL_777_B958, subset = nFeature_RNA >= 1500 & nFeature_RNA < 8000 & percent.mt < 10 & nCount_RNA > 2000 & nCount_RNA < 30000)
LCL_777_B958 <- NormalizeData(LCL_777_B958, normalization.method = "LogNormalize", scale.factor = 10000)
LCL_777_B958 <- FindVariableFeatures(LCL_777_B958, selection.method = "vst", nfeatures = 7000)
LCL_777_B958 <- ScaleData(LCL_777_B958, features = VariableFeatures(LCL_777_B958), verbose = TRUE)
LCL_777_B958 <- RunPCA(LCL_777_B958, features = VariableFeatures(object = LCL_777_B958))
LCL_777_B958 <- FindNeighbors(LCL_777_B958, dims = 1:20)
LCL_777_B958 <- FindClusters(LCL_777_B958, resolution = 0.1)
LCL_777_B958 <- RunUMAP(LCL_777_B958, dims = 1:20)
DimPlot(LCL_777_B958)

LCL_777_B958$orig.ident
get_expr_list(obj = LCL_777_B958, obj.ident = "orig.ident", output.name = paste0("4.", "LCL_777_B958", "_expr_list.rds"))





LCL_461_B958 <- Read10X(data.dir = "E:/Master_3rd_year/Projects/CellLines/PublicData/GSE158275_RAW_LCL/LCL_461_B958/")
LCL_461_B958 <- CreateSeuratObject(counts = LCL_461_B958, project = "LCL_461_B958", min.cells = 3, min.features = 200)
LCL_461_B958[["percent.mt"]] <- PercentageFeatureSet(LCL_461_B958, pattern = "^MT-")
VlnPlot(LCL_461_B958, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(LCL_461_B958, features = c("nFeature_RNA"), y.max = 2000)
VlnPlot(LCL_461_B958, features = c("nCount_RNA"), y.max = 8000)
VlnPlot(LCL_461_B958, features = c("percent.mt"), y.max = 10)

LCL_461_B958 <- subset(LCL_461_B958, subset = nFeature_RNA >= 1500 & nFeature_RNA < 8000 & percent.mt < 10 & nCount_RNA > 2000 & nCount_RNA < 30000)
LCL_461_B958 <- NormalizeData(LCL_461_B958, normalization.method = "LogNormalize", scale.factor = 10000)
LCL_461_B958 <- FindVariableFeatures(LCL_461_B958, selection.method = "vst", nfeatures = 7000)
LCL_461_B958 <- ScaleData(LCL_461_B958, features = VariableFeatures(LCL_461_B958), verbose = TRUE)
LCL_461_B958 <- RunPCA(LCL_461_B958, features = VariableFeatures(object = LCL_461_B958))
LCL_461_B958 <- FindNeighbors(LCL_461_B958, dims = 1:20)
LCL_461_B958 <- FindClusters(LCL_461_B958, resolution = 0.1)
LCL_461_B958 <- RunUMAP(LCL_461_B958, dims = 1:20)
DimPlot(LCL_461_B958)

LCL_461_B958$orig.ident
get_expr_list(obj = LCL_461_B958, obj.ident = "orig.ident", output.name = paste0("4.", "LCL_461_B958", "_expr_list.rds"))



LCL_777_B958_M81_AGGR <- Read10X(data.dir = "E:/Master_3rd_year/Projects/CellLines/PublicData/GSE158275_RAW_LCL/LCL_777_B958_M81_AGGR/")
LCL_777_B958_M81_AGGR <- CreateSeuratObject(counts = LCL_777_B958_M81_AGGR, project = "LCL_777_B958_M81_AGGR", min.cells = 3, min.features = 200)
LCL_777_B958_M81_AGGR[["percent.mt"]] <- PercentageFeatureSet(LCL_777_B958_M81_AGGR, pattern = "^MT-")
VlnPlot(LCL_777_B958_M81_AGGR, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(LCL_777_B958_M81_AGGR, features = c("nFeature_RNA"), y.max = 2000)
VlnPlot(LCL_777_B958_M81_AGGR, features = c("nCount_RNA"), y.max = 8000)
VlnPlot(LCL_777_B958_M81_AGGR, features = c("percent.mt"), y.max = 10)

LCL_777_B958_M81_AGGR <- subset(LCL_777_B958_M81_AGGR, subset = nFeature_RNA >= 1500 & nFeature_RNA < 8000 & percent.mt < 10 & nCount_RNA > 2000 & nCount_RNA < 30000)
LCL_777_B958_M81_AGGR <- NormalizeData(LCL_777_B958_M81_AGGR, normalization.method = "LogNormalize", scale.factor = 10000)
LCL_777_B958_M81_AGGR <- FindVariableFeatures(LCL_777_B958_M81_AGGR, selection.method = "vst", nfeatures = 7000)
LCL_777_B958_M81_AGGR <- ScaleData(LCL_777_B958_M81_AGGR, features = VariableFeatures(LCL_777_B958_M81_AGGR), verbose = TRUE)
LCL_777_B958_M81_AGGR <- RunPCA(LCL_777_B958_M81_AGGR, features = VariableFeatures(object = LCL_777_B958_M81_AGGR))
LCL_777_B958_M81_AGGR <- FindNeighbors(LCL_777_B958_M81_AGGR, dims = 1:20)
LCL_777_B958_M81_AGGR <- FindClusters(LCL_777_B958_M81_AGGR, resolution = 0.1)
LCL_777_B958_M81_AGGR <- RunUMAP(LCL_777_B958_M81_AGGR, dims = 1:20)
DimPlot(LCL_777_B958_M81_AGGR)

LCL_777_B958_M81_AGGR$orig.ident
get_expr_list(obj = LCL_777_B958_M81_AGGR, obj.ident = "orig.ident", output.name = paste0("4.", "LCL_777_B958_M81_AGGR", "_expr_list.rds"))



rm(LCL3)






#### Zhao - 2019 ####

Zhao <- read.table(paste0(public_dir, "GSM3689513_Mix4.expression_matrix.txt"))
Zhao[1:4, 1:4]

Zhao.seu <- CreateSeuratObject(Zhao, project = "Zhao_2019", min.cells = 3, min.features = 200)
Zhao.seu[["percent.mt"]] <- PercentageFeatureSet(Zhao.seu, pattern = "^MT-")

rb.genes <- rownames(Zhao.seu)[grep("^RP[SL]",rownames(Zhao.seu))]
Zhao.seu.count <- GetAssayData(object = Zhao.seu, slot = "counts")
percent.ribo <- Matrix::colSums(Zhao.seu.count[rb.genes,])/Matrix::colSums(Zhao.seu.count)*100
Zhao.seu <- AddMetaData(Zhao.seu, percent.ribo, col.name = "percent.ribo")
Zhao.seu_f <- subset(Zhao.seu, subset = (nFeature_RNA >= 1500 & nFeature_RNA <= 8000 & percent.mt <= 25 & nCount_RNA >= 2000 & nCount_RNA < 30000))
Zhao.seu_f
VlnPlot(Zhao.seu, c("nCount_RNA","nFeature_RNA","percent.mt"), group.by = "orig.ident", pt.size=0)

DefaultAssay(Zhao.seu_f) <- "RNA"
Zhao.seu_f <- NormalizeData(Zhao.seu_f, normalization.method = "LogNormalize", scale.factor = 10000)
Zhao.seu_f <- FindVariableFeatures(Zhao.seu_f, selection.method = "vst", nfeatures = 8000)
Zhao.seu_f <- ScaleData(Zhao.seu_f, features = VariableFeatures(Zhao.seu_f), verbose = TRUE)
Zhao.seu_f <- RunPCA(Zhao.seu_f, features = VariableFeatures(object = Zhao.seu_f))
Zhao.seu_f <- FindNeighbors(Zhao.seu_f, dims = 1:20)
Zhao.seu_f <- FindClusters(Zhao.seu_f, resolution = 0.1)
Zhao.seu_f <- RunUMAP(Zhao.seu_f, dims = 1:20)

DimPlot(Zhao.seu_f)




Zhao.3 <- read.table(paste0(public_dir, "GSM3689514_Mix3.expression_matrix.txt"))
Zhao.3[1:4, 1:4]
Zhao.seu3 <- CreateSeuratObject(Zhao.3, project = "Zhao_2019_mix3", min.cells = 3, min.features = 200)
Zhao.seu3
Zhao.seu3[["percent.mt"]] <- PercentageFeatureSet(Zhao.seu3, pattern = "^MT-")
rb.genes <- rownames(Zhao.seu3)[grep("^RP[SL]",rownames(Zhao.seu3))]
Zhao.seu3.count <- GetAssayData(object = Zhao.seu3, slot = "counts")
percent.ribo <- Matrix::colSums(Zhao.seu3.count[rb.genes,])/Matrix::colSums(Zhao.seu3.count)*100
Zhao.seu3 <- AddMetaData(Zhao.seu3, percent.ribo, col.name = "percent.ribo")
Zhao.seu3_f <- subset(Zhao.seu3, subset = (nFeature_RNA >= 1500 & nFeature_RNA <= 8000 & percent.mt <= 25 & nCount_RNA >= 2000 & nCount_RNA < 30000))
Zhao.seu3_f

DefaultAssay(Zhao.seu3_f) <- "RNA"
Zhao.seu3_f <- NormalizeData(Zhao.seu3_f, normalization.method = "LogNormalize", scale.factor = 10000)
Zhao.seu3_f <- FindVariableFeatures(Zhao.seu3_f, selection.method = "vst", nfeatures = 8000)
Zhao.seu3_f <- ScaleData(Zhao.seu3_f, features = VariableFeatures(Zhao.seu3_f), verbose = TRUE)
Zhao.seu3_f <- RunPCA(Zhao.seu3_f, features = VariableFeatures(object = Zhao.seu3_f))
Zhao.seu3_f <- FindNeighbors(Zhao.seu3_f, dims = 1:20)
Zhao.seu3_f <- FindClusters(Zhao.seu3_f, resolution = 0.1)
Zhao.seu3_f <- RunUMAP(Zhao.seu3_f, dims = 1:20)

DimPlot(Zhao.seu3_f)



####################################################################################-

#### 5.3 Shell ####

# input1: 4.CLL_expr_list.rds
# input2: Order of Cell Line (From 1 to 45)

args <- commandArgs(T)
if(length(args) != 2){
  print('Usage: NMF_auto.R 4.CLL_expr_list.rds Order')
  quit('Wrong')
}

library(NMF)

# input argument preprocess
expr_list <- readRDS(args[1])
i <- as.numeric(args[2])


print(paste0("####### Processing ",names(expr_list)[i]," ########"))
w <- NULL
h <- NULL
for(j in 5:15) {
  print(paste0("####### K factors at ",j," ########"))
  nmf_programs <- NMF::nmf(expr_list[[i]], rank=j, method="snmf/r", seed=123)
  saveRDS(nmf_programs, paste0("./nmf_programs_", strsplit(args[1],"\\.")[[1]][2], "_K", j, ".rds"))
  nmf_programs_scores <- list(w_basis=NMF::basis(nmf_programs), h_coef=t(NMF::coef(nmf_programs)))
  colnames(nmf_programs_scores$w_basis) <- paste0(names(expr_list)[i], "_", j, ".", 1:j)
  colnames(nmf_programs_scores$h_coef) <- paste0(names(expr_list)[i], "_", j, ".", 1:j)
  w <- cbind(w, nmf_programs_scores$w_basis)
  h <- cbind(h, nmf_programs_scores$h_coef)
}
saveRDS(w, paste0("./w_basis_", names(expr_list)[i], ".rds"))
saveRDS(h, paste0("./h_coef_", names(expr_list)[i], ".rds"))




#### 5.5 多线程后形成列表 ####

###### W & H
w_basis <- list() # nmf gene scores
h_coef <- list() # nmf cell scores

list.files(NMFout_dir, recursive = T) %>% length()
name.ls <- list.files(NMFout_dir, recursive = T, pattern = "w_basis")
CL.ls <- lapply(name.ls, function(x){strsplit(x, split = "/")}[[1]][1]) %>% unlist()
name.ls <- lapply(name.ls, function(x){strsplit(x, split = "/")}[[1]][2]) %>% unlist()
name.ls <- lapply(name.ls, function(x){strsplit(x, split = ".rds")}[[1]][1]) %>% unlist()
name.ls <- lapply(name.ls, function(x){strsplit(x, split = "_")}[[1]][3]) %>% unlist()
name.ls <- paste(CL.ls, name.ls, sep = "_")
head(name.ls); length(name.ls)

for (i in 1:length(name.ls)) {
  w <- readRDS(paste0(NMFout_dir, list.files(NMFout_dir, recursive = T, pattern = "w_basis")[i]))
  ll <- strsplit(colnames(w), split = "_")[[1]] %>% length()
  tt <- lapply(colnames(w), function(x){strsplit(x, split = "_")}[[1]][ll]) %>% unlist()
  colnames(w) <- paste0(name.ls[i], "_", tt)
  h <- readRDS(paste0(NMFout_dir, list.files(NMFout_dir, recursive = T, pattern = "h_coef")[i]))
  tt <- lapply(colnames(h), function(x){strsplit(x, split = "_")}[[1]][ll]) %>% unlist()
  colnames(h) <- paste0(name.ls[i], "_", tt)
  w_basis[[name.ls[i]]] <- w
  h_coef[[name.ls[i]]] <- h
}
names(w_basis) %>% order()


###### expr list

lapply(names(w_basis), function(x){strsplit(x, split = "_")}[[1]][1]) %>% unlist() %>% unique()
#### public data
fn.ls <- list.files(NMFin_dir, recursive = T)[c(1:11, 16:18)]
expr_list <- list()
for (i in 1:length(fn.ls)) {
  print(fn.ls[i])
  fn <- strsplit(fn.ls[i], split = "\\.")[[1]][2] # 取两个点之间的所属细胞系批次信息
  print(fn)
  test <- readRDS(paste0(NMFin_dir, fn.ls[i]))
  if (fn=="LCL") {
    names(test) <- paste0(strsplit(names(test), "_")[[1]], collapse = "")
  } else {
    names(test) <- lapply(names(test), function(x){strsplit(x, split = "_")}[[1]][1])
  }
  names(test) <- paste0(fn, "_", names(test))
  if ("CCLE_SCC25" %in% names(test)) {
    pooln <- strsplit(fn.ls[i], split = "\\.")[[1]][4]
    names(test)[names(test)=="CCLE_SCC25"] <- paste0("CCLE_SCC25", pooln)
  }
  expr_list <- append(expr_list, test)
}



#### CL
fn.ls <- list.files(NMFin_dir, recursive = T)[12:14]
lapply(fn.ls, function(x){strsplit(x, split = "\\.")}[[1]][2]) %>% unlist()
lapply(fn.ls, function(x){strsplit(x, split = "_")}[[1]][1]) %>% unlist()

for (i in 1:length(fn.ls)) {
  print(fn.ls[i])
  fn <- strsplit(fn.ls[i], split = "\\.")[[1]][2]
  fn <- strsplit(fn, split = "_")[[1]][1]
  print(fn)
  test <- readRDS(paste0(NMFin_dir, fn.ls[i]))
  names(test) <- paste0(fn, "_", names(test))
  expr_list <- append(expr_list, test)
}
names(expr_list)

expr_list <- expr_list[names(expr_list)!="CL1_MRC5"]
expr_list <- expr_list[names(expr_list)!="CL1_MCF10A"]
expr_list <- expr_list[names(expr_list) %>% order()]

identical(names(expr_list), names(w_basis))

for (i in names(expr_list)) {
  if (str_detect(i, "MCF10A")) {
    print(i)
  }
}

w_basis <- w_basis[names(w_basis)!="BC7_MCF10A"]
w_basis <- w_basis[names(w_basis)!="CL1_MRC5"]
w_basis <- w_basis[names(w_basis)!="CL1_MCF10A"]
h_coef <- h_coef[names(h_coef)!="BC7_MCF10A"]
h_coef <- h_coef[names(h_coef)!="CL1_MRC5"]
h_coef <- h_coef[names(h_coef)!="CL1_MCF10A"]
expr_list <- expr_list[names(expr_list)!="BC7_MCF10A"]

identical(names(expr_list), names(w_basis))
identical(names(expr_list), names(h_coef))

names(expr_list)[!(names(expr_list) %in% names(w_basis))]
length(expr_list)

rownames(w_basis$LCL_LCL461B958) <- lapply(rownames(w_basis$LCL_LCL461B958), function(x){strsplit(x, "GRCh38-")[[1]][2]}) %>% unlist
rownames(h_coef$LCL_LCL461B958) <- lapply(rownames(h_coef$LCL_LCL461B958), function(x){strsplit(x, "GRCh38-")[[1]][2]}) %>% unlist
rownames(expr_list$LCL_LCL461B958) <- lapply(rownames(expr_list$LCL_LCL461B958), function(x){strsplit(x, "GRCh38-")[[1]][2]}) %>% unlist


w_basis <- w_basis[-which(names(w_basis)=="CL3_A431")]
h_coef <- h_coef[-which(names(h_coef)=="CL3_A431")]
expr_list <- expr_list[-which(names(expr_list)=="CL3_A431")]
identical(names(expr_list), names(w_basis))
identical(names(expr_list), names(h_coef))
length(expr_list)



colnames(w_basis$BC32_AU565)





#### 5.6 NMF programs ####
## **(In this section, some codes were cited from the paper "Pan-cancer single-cell RNA-seq identifies recurring programs of cellular heterogeneity" and "Hallmarks of transcriptional intratumour heterogeneity across a thousand tumours")**

nmf_programs_sig <- lapply(w_basis, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:50]))
nmf_programs_sig$CCLE_2313287[1:4,1:4]

robust_nmf_programs <- function(nmf_programs, intra_min = 35, intra_max = 10, inter_filter=T, inter_min = 10) {
  intra_intersect <- lapply(nmf_programs, function(z) apply(z, 2, function(x) apply(z, 2, function(y) length(intersect(x,y)))))
  intra_intersect_max <- lapply(intra_intersect, function(x) apply(x, 2, function(y) sort(y, decreasing = T)[2]))
  nmf_sel <- lapply(names(nmf_programs), function(x) nmf_programs[[x]][,intra_intersect_max[[x]]>=intra_min])
  names(nmf_sel) <- names(nmf_programs)
  nmf_sel_unlist <- do.call(cbind, nmf_sel)
  inter_intersect <- apply(nmf_sel_unlist , 2, function(x) apply(nmf_sel_unlist , 2, function(y) length(intersect(x,y))))
  
  final_filter <- NULL 
  for(i in names(nmf_sel)) {
    a <- inter_intersect[grep(i, colnames(inter_intersect), invert = T),grep(i, colnames(inter_intersect))]
    b <- sort(apply(a, 2, max), decreasing = T)
    if(inter_filter==T) b <- b[b>=inter_min]
    if(length(b) > 1) {
      c <- names(b[1])
      for (y in 2:length(b)) { 
        if(max(inter_intersect[c,names(b[y])]) <= intra_max) c <- c(c,names(b[y]))
      }
      final_filter <- c(final_filter, c)
    } else {
      final_filter <- c(final_filter, names(b))
    }
  }
  return(final_filter)
}
nmf_filter <- robust_nmf_programs(nmf_programs_sig, intra_min = 35, intra_max = 10, inter_filter=T, inter_min = 10)

nmf_programs_sig_f <- lapply(nmf_programs_sig, function(x) x[, is.element(colnames(x), nmf_filter),drop=F])
nmf_programs_sig_f  <- do.call(cbind, nmf_programs_sig_f)
nmf_intersect <- apply(nmf_programs_sig_f , 2, function(x) apply(nmf_programs_sig_f , 2, function(y) length(intersect(x,y))))
dim(nmf_intersect)
nmf_intersect[1:20, 1:4]
nmf_intersect_hc <- hclust(as.dist(50-nmf_intersect), method="average")
nmf_intersect_hc <- reorder(as.dendrogram(nmf_intersect_hc), colMeans(nmf_intersect))
nmf_intersect <- nmf_intersect[order.dendrogram(nmf_intersect_hc), order.dendrogram(nmf_intersect_hc)]

nmf_intersect_meltI <- reshape2::melt(nmf_intersect)

library(RColorBrewer)
library(viridis)
custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

p2.1 <- ggplot(data = nmf_intersect_meltI, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) +
  geom_tile() +
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(),
         panel.border = element_rect(fill=F),
         panel.background = element_blank(),
         axis.line = element_blank(),
         axis.text = element_blank(),
         axis.title = element_blank(),
         legend.title = element_text(size=11),
         legend.text = element_text(size = 10),
         legend.text.align = 0.5,
         legend.justification = "bottom") +
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))



#### 5.7 corr ####
## **(In this section, some codes were cited from the paper "Pan-cancer single-cell RNA-seq identifies recurring programs of cellular heterogeneity" and "Hallmarks of transcriptional intratumour heterogeneity across a thousand tumours")**

nmf_programs_cells <- h_coef
complexity_expr <- lapply(expr_list, function(x) apply(x, 2, function(y) length(which(y != 0))))   
nmf_corr_complexity_expr <- data.frame(matrix(ncol = 1, nrow = ncol(nmf_intersect)), row.names=colnames(nmf_intersect))
colnames(nmf_corr_complexity_expr) <- "corr"
for(i in rownames(nmf_corr_complexity_expr)) {
  a <- paste(strsplit(i, "_")[[1]][1:2], collapse = "_")
  b <- nmf_programs_cells[[a]][,i]
  c <- complexity_expr[[a]]
  nmf_corr_complexity_expr[i,1] <- cor(b,c)
}
head(nmf_corr_complexity_expr)
p2.2 <- ggplot(nmf_corr_complexity_expr, aes(y=corr, x = 1:nrow(nmf_corr_complexity_expr))) +
  geom_smooth(span=0.05, se = FALSE, color="gray36", size= 0.8, method = "loess") + 
  geom_point(alpha=0.3, size = 1.5) + 
  theme(axis.ticks = element_blank(), axis.line = element_blank(), panel.background = element_rect(fill = "gray97"), panel.grid = element_blank(), axis.text.x = element_blank(), axis.title = element_text(size = 12), axis.text.y = element_text(size = 10), plot.margin=unit(c(1,1,-0.6,1), "cm") ) + 
  labs(y="Cell complexity\ncorrelation", x="") + 
  scale_y_continuous(expand = c(0,0), limits = c(-0.9,0.9), breaks = seq(-0.8, 0.8, 0.4)) +
  scale_x_continuous(expand = c(0,0))


#### 5.8 filtering ####

corr_cutoff <- 0.5
nmf_corr_complexity_expr_f <- subset(nmf_corr_complexity_expr, corr < corr_cutoff & corr > -corr_cutoff)

nmf_programs_sig_ff <- nmf_programs_sig_f[,which(colnames(nmf_programs_sig_f) %in% rownames(nmf_corr_complexity_expr_f))]

nmf_intersect_f <- apply(nmf_programs_sig_ff , 2, function(x) apply(nmf_programs_sig_ff , 2, function(y) length(intersect(x,y))))
nmf_intersect_hc <- hclust(as.dist(50-nmf_intersect_f), method="average")
nmf_intersect_hc <- reorder(as.dendrogram(nmf_intersect_hc), colMeans(nmf_intersect_f))
nmf_intersect_f <- nmf_intersect_f[order.dendrogram(nmf_intersect_hc), order.dendrogram(nmf_intersect_hc)]
nmf_intersect_f[1:4,1:4]

nmf_intersect <- nmf_intersect_f



#### 5.9 clustering ####
## **(In this section, some codes were cited from the paper "Pan-cancer single-cell RNA-seq identifies recurring programs of cellular heterogeneity" and "Hallmarks of transcriptional intratumour heterogeneity across a thousand tumours")**

nmf_programs <- nmf_programs_sig_ff
Genes_nmf_w_basis <- w_basis

Min_intersect_initial <- 10
Min_intersect_cluster <- 10
Min_group_size        <- 5

Sorted_intersection <- sort(apply(nmf_intersect , 2, function(x) (length(which(x>=Min_intersect_initial))-1)  ) , decreasing = TRUE)
Sorted_intersection %>% head
identical(colnames(nmf_programs_sig_f) %>% sort, names(Sorted_intersection) %>% sort)

Cluster_list              <- list()
MP_list                   <- list()
k                         <- 1
Curr_cluster              <- c()

nmf_intersect_original    <- nmf_intersect

while (Sorted_intersection[1]>Min_group_size) {  
  
  Curr_cluster <- c(Curr_cluster , names(Sorted_intersection[1]))
  Genes_MP                    <- nmf_programs[,names(Sorted_intersection[1])]
  nmf_programs                <- nmf_programs[,-match(names(Sorted_intersection[1]) , colnames(nmf_programs))]
  nmf_programs %>% colnames %>% length
  Intersection_with_Genes_MP  <- sort(apply(nmf_programs, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE)
  Intersection_with_Genes_MP %>% head
  NMF_history                 <- Genes_MP
  while ( Intersection_with_Genes_MP[1] >= Min_intersect_cluster) {
    
    Curr_cluster  <- c(Curr_cluster , names(Intersection_with_Genes_MP)[1])
    Curr_cluster
    Genes_MP_temp   <- sort(table(c(NMF_history , nmf_programs[,names(Intersection_with_Genes_MP)[1]])), decreasing = TRUE)
    Genes_MP_temp
    Genes_at_border <- Genes_MP_temp[which(Genes_MP_temp == Genes_MP_temp[50])]
    Genes_at_border
    
    if (length(Genes_at_border)>1){
      Genes_curr_NMF_score <- c()
      for (i in Curr_cluster) {
        curr_study           <- paste( strsplit(i , "[_]")[[1]][1:2]   , collapse = "_"  )
        curr_study
        Q                    <- Genes_nmf_w_basis[[curr_study]][ match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]])))[!is.na(match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]]))))]   ,i] 
        Q
        names(Q)             <- names(Genes_at_border[!is.na(match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]]))))])
        Genes_curr_NMF_score <- c(Genes_curr_NMF_score,  Q )
      }
      Genes_curr_NMF_score_sort <- sort(Genes_curr_NMF_score , decreasing = TRUE)
      Genes_curr_NMF_score_sort <- Genes_curr_NMF_score_sort[unique(names(Genes_curr_NMF_score_sort))]
      Genes_MP_temp             <- c(names(Genes_MP_temp[which(Genes_MP_temp > Genes_MP_temp[50])]) , names(Genes_curr_NMF_score_sort))
      
    } else {
      Genes_MP_temp <- names(Genes_MP_temp)[1:50] 
    }
    NMF_history     <- c(NMF_history , nmf_programs[,names(Intersection_with_Genes_MP)[1]]) 
    Genes_MP        <- Genes_MP_temp[1:50]
    nmf_programs    <- nmf_programs[,-match(names(Intersection_with_Genes_MP)[1] , colnames(nmf_programs))] 
    Intersection_with_Genes_MP <- sort(apply(nmf_programs, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE)
    
  }
  Cluster_list[[paste0("Cluster_", k)]] <- Curr_cluster
  MP_list[[paste0("MP_", k)]]           <- Genes_MP
  k <- k + 1
  nmf_intersect             <- nmf_intersect[-match(Curr_cluster,rownames(nmf_intersect) ) , -match(Curr_cluster,colnames(nmf_intersect) ) ]
  dim(nmf_intersect)
  Sorted_intersection       <-  sort(apply(nmf_intersect , 2, function(x) (length(which(x>=Min_intersect_initial))-1)  ) , decreasing = TRUE)
  Sorted_intersection
  Curr_cluster <- c()
  print(dim(nmf_intersect)[2])
}

nmf_intersect <- nmf_intersect_original



#### 5.10 plot ####

# Custom color palette
library(RColorBrewer)
library(viridis)
custom_magma <- c(colorRampPalette(c("white", rev(magma(323, begin = 0.15))[1]))(10), rev(magma(323, begin = 0.18)))

Cluster_list2 <- Cluster_list[order(sapply(Cluster_list, function(x) length(x)), decreasing = T)]
names(Cluster_list2) <- paste0("Cluster_", 1:length(Cluster_list2))
Cluster_list2

MP_list2 <- MP_list
MP_list2 <- MP_list2[order(sapply(Cluster_list, function(x) length(x)), decreasing = T)]
names(MP_list2) <- paste0("MP_", 1:length(MP_list2))
MP_list2

MP_list.id <- MP_list2
for (i in 1:length(MP_list2)) {
  MP_list.id[[i]] <- bitr(MP_list2[[i]], fromType = 'SYMBOL', toType = c('ENTREZID'), OrgDb = org.Hs.eg.db)$ENTREZID
}
MP_list.id


#### 5.10.1 arrange ####

#### RHP
RHP.data <- readxl::read_xlsx("./IF38.33.NG.41588_2020_726_MOESM3_ESM.xlsx", sheet = "Table S4", skip=3)
RHP.df <- tibble(RHP="", gene="")
for (i in 1:ncol(RHP.data)) {
  test2 <- na.omit(c(RHP.data[,i])[[1]])
  test <- tibble(RHP=rep(colnames(RHP.data)[i], length(test2)), gene=test2)
  RHP.df <- rbind(RHP.df, test)
}
RHP.df <- RHP.df[-1,]
RHP.df

Barkley.module <- readxl::read_xlsx(path = "./IF41.307.NatureGenetic.Barkley2022_41588_2022_1141_MOESM6_ESM.xlsx", sheet="TableS3")
colnames(Barkley.module) <- Barkley.module[2,]
Barkley.module <- Barkley.module[-c(1,2),]
Barkley.df <- tibble(Barkley="", gene="")
for (i in 1:ncol(Barkley.module)) {
  test <- tibble(Barkley=rep(colnames(Barkley.module)[i], length(c(Barkley.module[,i])[[1]])), gene=c(Barkley.module[,i])[[1]])
  Barkley.df <- rbind(Barkley.df, test)
}
Barkley.df <- Barkley.df[-1,]
Barkley.df

Nat.MP <- readxl::read_xlsx("./IF69.504.Nature.Gavish2023_41586_2023_6130_MOESM6_ESM.xlsx")
Nat.MP.df <- tibble(Nat.MP="", gene="")
for (i in 1:ncol(Nat.MP)) {
  test2 <- c(Nat.MP[,i])[[1]]
  test <- tibble(Nat.MP=rep(colnames(Nat.MP)[i], length(test2)), gene=test2)
  Nat.MP.df <- rbind(Nat.MP.df, test)
}
Nat.MP.df <- Nat.MP.df[-1,]
Nat.MP.df

library(msigdbr)
Hs.H <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, human_gene_symbol)


NC.rp <- readxl::read_xlsx("./IF16.6.NC.Zhu2023_41467_2023_43991_MOESM5_ESM.xlsx", skip = 2)
NC.rp.df <- tibble(NC.rp="", gene="")
for (i in 1:nrow(NC.rp)) {
  test2 <- lapply(NC.rp$`Frequent  genes`[i], function(x){strsplit(x, " ")}) %>% unlist
  test <- tibble(NC.rp=rep(NC.rp$`NMF based program`[i], length(test2)), gene=test2)
  NC.rp.df <- rbind(NC.rp.df, test)
}
NC.rp.df <- NC.rp.df[-1,]
NC.rp.df



enrich.df <- tibble(MP="", RHP="", Barkley="", Nat3CA="", HM="", NC="")

for (i in names(MP_list2)) {
  
  tryCatch(
    {go <- enricher(MP_list2[[i]], TERM2GENE = RHP.df)},
    error = function(e){
      message(paste0(i, ' HAS NO RESULTS'))
      print(e)
      return(NULL)
    },
    warning = function(w){
      message('A warning')
      print(w)
      return(NULL)
    }
  )
  if (is.null(go)) {
    RHP.go <- ""
  } else if (dim(go)[1]!=0) {
    RHP.go <- go@result$Description[1]
  } else {
    RHP.go <- ""
  }
  
  ####
  tryCatch(
    {go <- enricher(MP_list2[[i]], TERM2GENE = Barkley.df)},
    error = function(e){
      message(paste0(i, ' HAS NO RESULTS'))
      print(e)
      return(NULL)
    },
    warning = function(w){
      message('A warning')
      print(w)
      return(NULL)
    }
  )
  if (is.null(go)) {
    Barkley.go <- ""
  } else if (dim(go)[1]!=0) {
    Barkley.go <- go@result$Description[1]
  } else {
    Barkley.go <- ""
  }
  
  ####
  tryCatch(
    {go <- enricher(MP_list2[[i]], TERM2GENE = Nat.MP.df)},
    error = function(e){
      message(paste0(i, ' HAS NO RESULTS'))
      print(e)
      return(NULL)
    },
    warning = function(w){
      message('A warning')
      print(w)
      return(NULL)
    }
  )
  if (is.null(go)) {
    Nat.MP.go <- ""
  } else if (dim(go)[1]!=0) {
    Nat.MP.go <- go@result$Description[1]
  } else {
    Nat.MP.go <- ""
  }
  
  ####
  tryCatch(
    {go <- enricher(MP_list2[[i]], TERM2GENE = Hs.H)},
    error = function(e){
      message(paste0(i, ' HAS NO RESULTS'))
      print(e)
      return(NULL)
    },
    warning = function(w){
      message('A warning')
      print(w)
      return(NULL)
    }
  )
  if (is.null(go)) {
    HM.go <- ""
  } else if (dim(go)[1]!=0) {
    HM.go <- go@result$Description[1]
  } else {
    HM.go <- ""
  }
  
  ####
  tryCatch(
    {go <- enricher(MP_list2[[i]], TERM2GENE = NC.rp.df)},
    error = function(e){
      message(paste0(i, ' HAS NO RESULTS'))
      print(e)
      return(NULL)
    },
    warning = function(w){
      message('A warning')
      print(w)
      return(NULL)
    }
  )
  if (is.null(go)) {
    NC.go <- ""
  } else if (dim(go)[1]!=0) {
    NC.go <- go@result$Description[1]
  } else {
    NC.go <- ""
  }
  
  test.df <- tibble(MP=i, RHP=RHP.go, Barkley=Barkley.go, Nat3CA=Nat.MP.go, HM=HM.go, NC=NC.go)
  enrich.df <- rbind(enrich.df, test.df)
}
enrich.df <- enrich.df[-1,]
enrich.df
as.data.frame(enrich.df)
as.data.frame(enrich.df[,c("MP", "RHP", "Nat3CA", "NC")])




idx <- c(1, 10, 20, # G2M
         2:4, # G1S, Proteasomal degradation, Stress
         18, 5, 11, 17, 8, # EMT 1 2 2 3 4
  6, # Interferon / MHC-2
  7, # Estrogen response
  19, # Epithelial Senescence
  9, # Secreted
  12, 13, # Chromatin
  14:16, # Hypoxia/chol, Translation initiation, Skin-pigmentation
  21 # p53-Dependent Senescence
  )
length(unique(idx))==length(Cluster_list2)


Cluster_list3 <- Cluster_list2[idx]
names(Cluster_list3) <- paste0("Cluster_", 1:length(Cluster_list3))

MP_list3 <- MP_list2[idx]
names(MP_list3) <- paste0("MP_", 1:length(MP_list3))

MP_list.id <- MP_list3
for (i in 1:length(MP_list3)) {
  MP_list.id[[i]] <- bitr(MP_list3[[i]], fromType = 'SYMBOL', toType = c('ENTREZID'), OrgDb = org.Hs.eg.db)$ENTREZID
}
MP_list.id


enrich.df <- tibble(MP="", RHP="", Barkley="", Nat3CA="", HM="", NC="")

for (i in names(MP_list3)) {
  
  tryCatch(
    {go <- enricher(MP_list3[[i]], TERM2GENE = RHP.df)},
    error = function(e){
      message(paste0(i, ' HAS NO RESULTS'))
      print(e)
      return(NULL)
    },
    warning = function(w){
      message('A warning')
      print(w)
      return(NULL)
    }
  )
  if (is.null(go)) {
    RHP.go <- ""
  } else if (dim(go)[1]!=0) {
    RHP.go <- go@result$Description[1]
  } else {
    RHP.go <- ""
  }
  
  tryCatch(
    {go <- enricher(MP_list3[[i]], TERM2GENE = Barkley.df)},
    error = function(e){
      message(paste0(i, ' HAS NO RESULTS'))
      print(e)
      return(NULL)
    },
    warning = function(w){
      message('A warning')
      print(w)
      return(NULL)
    }
  )
  if (is.null(go)) {
    Barkley.go <- ""
  } else if (dim(go)[1]!=0) {
    Barkley.go <- go@result$Description[1]
  } else {
    Barkley.go <- ""
  }
  
  tryCatch(
    {go <- enricher(MP_list3[[i]], TERM2GENE = Nat.MP.df)},
    error = function(e){
      message(paste0(i, ' HAS NO RESULTS'))
      print(e)
      return(NULL)
    },
    warning = function(w){
      message('A warning')
      print(w)
      return(NULL)
    }
  )
  if (is.null(go)) {
    Nat.MP.go <- ""
  } else if (dim(go)[1]!=0) {
    Nat.MP.go <- go@result$Description[1]
  } else {
    Nat.MP.go <- ""
  }
  
  tryCatch(
    {go <- enricher(MP_list3[[i]], TERM2GENE = Hs.H)},
    error = function(e){
      message(paste0(i, ' HAS NO RESULTS'))
      print(e)
      return(NULL)
    },
    warning = function(w){
      message('A warning')
      print(w)
      return(NULL)
    }
  )
  if (is.null(go)) {
    HM.go <- ""
  } else if (dim(go)[1]!=0) {
    HM.go <- go@result$Description[1]
  } else {
    HM.go <- ""
  }
  
  ####
  tryCatch(
    {go <- enricher(MP_list3[[i]], TERM2GENE = NC.rp.df)},
    error = function(e){
      message(paste0(i, ' HAS NO RESULTS'))
      print(e)
      return(NULL)
    },
    warning = function(w){
      message('A warning')
      print(w)
      return(NULL)
    }
  )
  if (is.null(go)) {
    NC.go <- ""
  } else if (dim(go)[1]!=0) {
    NC.go <- go@result$Description[1]
  } else {
    NC.go <- ""
  }
  
  test.df <- tibble(MP=i, RHP=RHP.go, Barkley=Barkley.go, Nat3CA=Nat.MP.go, HM=HM.go, NC=NC.go)
  enrich.df <- rbind(enrich.df, test.df)
}
enrich.df <- enrich.df[-1,]
enrich.df
as.data.frame(enrich.df)
as.data.frame(enrich.df[,c("MP", "RHP", "Nat3CA", "NC")])




inds_sorted <- c()
for (j in 1:length(Cluster_list3)){
  inds_sorted <- c(inds_sorted , match(Cluster_list3[[j]] , colnames(nmf_intersect_original)))
}
inds_new <- c(inds_sorted,   which(is.na( match(1:dim(nmf_intersect_original)[2],inds_sorted))))


nmf_intersect_meltI_NEW <- reshape2::melt(nmf_intersect_original[rev(inds_sorted), rev(inds_sorted)])
nmf_intersect_meltI_NEW$group <- "Cluster_Other"
for (i in names(Cluster_list3)) {
  nmf_intersect_meltI_NEW$group[nmf_intersect_meltI_NEW$Var1 %in% Cluster_list3[[i]]] <- i
}
nmf_intersect_meltI_NEW$group <- factor(nmf_intersect_meltI_NEW$group, levels = mixedsort(unique(nmf_intersect_meltI_NEW$group)))
df2 <- tibble(y=nmf_intersect_meltI_NEW$Var1) %>% as.data.frame()
df2$x <- 1
df2$group <- nmf_intersect_meltI_NEW$group
p2 <- ggplot(df2, aes(x=x,y=y))+
  geom_tile(aes(fill=group))+
  scale_x_continuous(expand = c(0,0))+
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "left",
        legend.title = element_blank()) +
  scale_fill_manual(values = c(fig.colors, dark.fig.colors))
p1 <- ggplot(data = nmf_intersect_meltI_NEW, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +                                
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))
p3 <- p1 %>% insert_left(p2,width = 0.05)




#### Include Others
nmf_intersect_meltI_NEW <- reshape2::melt(nmf_intersect_original[rev(inds_new), rev(inds_new)])
nmf_intersect_meltI_NEW$group <- "Cluster_Other"
for (i in names(Cluster_list3)) {
  nmf_intersect_meltI_NEW$group[nmf_intersect_meltI_NEW$Var1 %in% Cluster_list3[[i]]] <- i
}
nmf_intersect_meltI_NEW$group <- factor(nmf_intersect_meltI_NEW$group, levels = mixedsort(unique(nmf_intersect_meltI_NEW$group)))
df2 <- tibble(y=nmf_intersect_meltI_NEW$Var1) %>% as.data.frame()
df2$x <- 1
df2$group <- nmf_intersect_meltI_NEW$group
p2 <- ggplot(df2, aes(x=x,y=y))+
  geom_tile(aes(fill=group))+
  scale_x_continuous(expand = c(0,0))+
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "left",
        legend.title = element_blank()) +
  scale_fill_manual(values = c(fig.colors[1:21], "grey"))
p1 <- ggplot(data = nmf_intersect_meltI_NEW, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +                                
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))
p3 <- p1 %>% insert_left(p2,width = 0.05)



#### 5.11 NMF programs - Enrichment ####
Cluster_list3
MP_list3

#### 5.11.1 GO BP ####
GOBP <- compareCluster(MP_list3, fun = "enrichGO", ont = "BP",OrgDb=org.Hs.eg.db, keyType = 'SYMBOL')
pGOBP <- dotplot_choose3(GOBP, choose = 'top', term = 4, multi.thd = T) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 70)) +
  labs(title = "GO BP")
pGOBP



#### 5.11.2 GO BP ####

for (MP.idx in 1:length(MP_list3)) {
  tryCatch(
    {GOBP <- clusterProfiler::compareCluster(geneClusters = split(MP_list3[[MP.idx]], names(MP_list3)[MP.idx]), 
                                             fun = "enrichGO", ont = "BP", OrgDb = org.Hs.eg.db, 
                                             keyType = 'SYMBOL')},
    error = function(e){
      message(paste0(names(MP_list3)[MP.idx], ' HAS NO RESULTS FOR GOBP'))
      print(e)
      return(NULL)
    },
    warning = function(w){
      message('A warning')
      print(w)
      return(NULL)
    }
  )
  if (is.null(GOBP)) {
    next
  }
  dim(GOBP)
  GOBP <- clusterProfiler::simplify(GOBP)
  dim(GOBP)
  test3 <- GOBP@compareClusterResult
  for (i in 1:nrow(test3)) {
    test3$Gene_ratio[i] <- as.numeric(strsplit(test3$GeneRatio[i],"/")[[1]][1])/as.numeric(strsplit(test3$GeneRatio[i],"/")[[1]][2])
  }
  test3 <- test3 %>% top_n(n=10,wt=-p.adjust)
  ggplot(test3,
         aes(x=-log(p.adjust),
             y=reorder(Description,-log(p.adjust)),fill=Gene_ratio)) + 
    geom_bar(stat="identity", width=.4) + 
    theme1() + 
    scale_fill_gradientn(colours = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)) + 
    labs(y="", title=names(MP_list3)[MP.idx]) + 
    scale_y_discrete(labels = function(x) str_wrap(x, width = 50) ) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(size = 12))
  ggsave(paste0("6.BarPlot_MP_", MP.idx, "_GOBP.png"),width = 7, height = 4)
  ggsave(paste0("6.BarPlot_MP_", MP.idx, "_GOBP.pdf"),width = 7, height = 4)
}



#### heatmap ####

hyper.hm.mtx <- tibble(MetaProgram="", RHP=0, Barkley=0, Nat3CA=0, HM=0, NC=0)
for (i in 1:nrow(enrich.df)) {
  test <- tibble(MetaProgram=names(MP_list3)[i], RHP=NA, Barkley=NA, Nat3CA=NA, HM=NA, NC=NA)
  
  ## RHP
  go <- enricher(MP_list3[[names(MP_list3)[i]]], TERM2GENE = RHP.df)
  if (enrich.df$RHP[i]!="") {
    test$RHP <- go@result$p.adjust[which(go@result$ID==enrich.df$RHP[i])]
  }
  ## Barkley
  go <- enricher(MP_list3[[names(MP_list3)[i]]], TERM2GENE = Barkley.df)
  if (enrich.df$Barkley[i]!="") {
    test$Barkley <- go@result$p.adjust[which(go@result$ID==enrich.df$Barkley[i])]
  }
  ## Nat3CA
  go <- enricher(MP_list3[[names(MP_list3)[i]]], TERM2GENE = Nat.MP.df)
  if (enrich.df$Nat3CA[i]!="") {
    test$Nat3CA <- go@result$p.adjust[which(go@result$ID==enrich.df$Nat3CA[i])]
  }
  ## HM
  go <- enricher(MP_list3[[names(MP_list3)[i]]], TERM2GENE = Hs.H)
  if (enrich.df$HM[i]!="") {
    test$HM <- go@result$p.adjust[which(go@result$ID==enrich.df$HM[i])]
  }
  ## NC
  go <- enricher(MP_list3[[names(MP_list3)[i]]], TERM2GENE = NC.rp.df)
  if (enrich.df$NC[i]!="") {
    test$NC <- go@result$p.adjust[which(go@result$ID==enrich.df$NC[i])]
  }
  
  hyper.hm.mtx <- rbind(hyper.hm.mtx, test)
}
hyper.hm.mtx <- hyper.hm.mtx[-1,]
hyper.hm.mtx

my_color = colorRampPalette(c(season.colors3[16], season.colors3[9]))(50)


test <- hyper.hm.mtx[,2:ncol(hyper.hm.mtx)]
test <- as.data.frame(test)
rownames(test) <- paste0("MP_",1:21)
# colnames(test)

anno.df <- tibble(`Meta Program`=paste0("MP_",1:21)) %>% as.data.frame()
rownames(anno.df) <- rownames(test)

anno.color <- list(`Meta Program` = c(MP_1=fig.colors[1],
                                      MP_2=fig.colors[2],
                                      MP_3=fig.colors[3],
                                      MP_4=fig.colors[4],
                                      MP_5=fig.colors[5],
                                      MP_6=fig.colors[6],
                                      MP_7=fig.colors[7],
                                      MP_8=fig.colors[8],
                                      MP_9=fig.colors[9],
                                      MP_10=fig.colors[10],
                                      MP_11=fig.colors[11],
                                      MP_12=fig.colors[12],
                                      MP_13=fig.colors[13],
                                      MP_14=fig.colors[14],
                                      MP_15=fig.colors[15],
                                      MP_16=fig.colors[16],
                                      MP_17=fig.colors[17],
                                      MP_18=fig.colors[18],
                                      MP_19=fig.colors[19],
                                      MP_20=fig.colors[20],
                                      MP_21=dark.fig.colors[1]))
pheatmap::pheatmap(test[1:5], 
                   na_col = "light grey", 
                   cluster_rows = F, 
                   cluster_cols = F,
                   color = rev(my_color), 
                   angle_col = 45,
                   border_color = "white", 
                   annotation_row = anno.df, 
                   show_rownames = F,
                   gaps_row = c(1:21), 
                   gaps_col = 1:5, annotation_legend = F, annotation_names_row = F,
                   annotation_colors = anno.color)
dev.off()
pheatmap::pheatmap(test[1:5], 
                   na_col = "light grey", 
                   cluster_rows = F, 
                   cluster_cols = F,
                   color = rev(my_color), 
                   angle_col = 45,
                   border_color = "white", 
                   annotation_row = anno.df, 
                   show_rownames = F,
                   gaps_row = c(1:21), 
                   gaps_col = 1:5, annotation_legend = F, annotation_names_row = F,
                   annotation_colors = anno.color)
dev.off()






#### 5.11.4 GoPlot ####
for (i in 20:length(MP_list.id)) {
  ego <- enrichGO(gene          = MP_list.id[[i]],
                  # universe      = names(MP_list.id)[i],
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  p1 <- goplot(ego, showCategory = 10, layout = "sugiyama") + labs(title = names(MP_list.id)[i])
  ggsave(paste0("./GoPlot/5.GoPlot_", names(MP_list.id)[i], ".png"), width = 15, height = 10)
  ggsave(paste0("./GoPlot/5.GoPlot_", names(MP_list.id)[i], ".pdf"), width = 15, height = 10)
}






#### 5.11.5 KEGG ####

KEGG.res <- compareCluster(geneCluster = MP_list.id, fun = "enrichKEGG", organism = 'hsa', keyType = 'kegg')
pKEGG.res <- dotplot_choose3(KEGG.res,choose = 'top', term = 4, multi.thd = T)+theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(title = "KEGG")


#### 5.11.6 KEGG ####

for (MP.idx in 1:length(MP_list.id)) {
  tryCatch(
    {KEGG.res <- clusterProfiler::compareCluster(geneClusters = split(MP_list.id[[MP.idx]], names(MP_list.id)[MP.idx]), 
                                             fun = "enrichKEGG")},
    error = function(e){
      message(paste0(names(MP_list.id)[MP.idx], ' HAS NO RESULTS FOR KEGG'))
      print(e)
      return(NULL)
    },
    warning = function(w){
      message('A warning')
      print(w)
      return(NULL)
    }
  )
  if (is.null(KEGG.res)) {
    next
  }
  dim(KEGG.res)
  dim(KEGG.res)
  test3 <- KEGG.res@compareClusterResult
  for (i in 1:nrow(test3)) {
    test3$Gene_ratio[i] <- as.numeric(strsplit(test3$GeneRatio[i],"/")[[1]][1])/as.numeric(strsplit(test3$GeneRatio[i],"/")[[1]][2])
  }
  test3 <- test3 %>% top_n(n=10,wt=-p.adjust)
  ggplot(test3,
         aes(x=-log(p.adjust),
             y=reorder(Description,-log(p.adjust)),fill=Gene_ratio)) + 
    geom_bar(stat="identity", width=.4) + 
    theme1() + 
    scale_fill_gradientn(colours = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)) + 
    labs(y="", title=names(MP_list.id)[MP.idx]) + 
    scale_y_discrete(labels = function(x) str_wrap(x, width = 50) ) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(size = 12))
}




#### 5.11.7 GenClip3 ####

for (i in 1:length(MP_list3)) {
  write.table(x = MP_list3[[i]], file = paste0("6.ToGenClip3_", names(MP_list3)[i], ".txt"), row.names = F, col.names = F)
}


#### 5.11.6 ReactomePA ####
library(ReactomePA)
reactome.res <- compareCluster(geneCluster = MP_list.id, fun = "enrichPathway", organism = 'human')
preactome.res <- dotplot_choose3(reactome.res, choose = 'top', term = 4, multi.thd = T) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Reactome")
preactome.res


#### 5.11.7 DO ####
library(DOSE)
DO.res <- compareCluster(geneCluster = MP_list.id, fun = "enrichDO")
# pDO.res <- dotplot_choose(DO.res, choose = 'top', term = 4)+theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   labs(title = "DO MP")
pDO.res <- dotplot_choose3(DO.res, choose = 'top', term = 4, multi.thd = T)+theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(title = "DOSE")
pDO.res




#### 5.12 plotting ####

# Style1
nmf_intersect_meltI_NEW <- reshape2::melt(nmf_intersect_original[rev(inds_sorted), rev(inds_sorted)])
nmf_intersect_meltI_NEW$group <- "Cluster_Other"
for (i in names(Cluster_list3)) {
  nmf_intersect_meltI_NEW$group[nmf_intersect_meltI_NEW$Var1 %in% Cluster_list3[[i]]] <- i
}
nmf_intersect_meltI_NEW$group <- factor(nmf_intersect_meltI_NEW$group, levels = mixedsort(unique(nmf_intersect_meltI_NEW$group)))
df2 <- tibble(y=nmf_intersect_meltI_NEW$Var1) %>% as.data.frame()
df2$x <- 1
df2$group <- nmf_intersect_meltI_NEW$group

df2$Program <- NA
unique(df2$group)
df2$Program[df2$group=="Cluster_1"] <- "CellCycle_G2M_1"
df2$Program[df2$group=="Cluster_2"] <- "CellCycle_G2M_2"
df2$Program[df2$group=="Cluster_3"] <- "CellCycle_G2M_3"
df2$Program[df2$group=="Cluster_4"] <- "CellCycle_G1S"
df2$Program[df2$group=="Cluster_5"] <- "Proteasomal_Degradation"
df2$Program[df2$group=="Cluster_6"] <- "Stress"
df2$Program[df2$group=="Cluster_7"] <- "EMT-I"
df2$Program[df2$group=="Cluster_8"] <- "EMT-II_1"
df2$Program[df2$group=="Cluster_9"] <- "EMT-II_2"
df2$Program[df2$group=="Cluster_10"] <- "EMT-III"
df2$Program[df2$group=="Cluster_11"] <- "EMT-IV"
df2$Program[df2$group=="Cluster_12"] <- "Interferon/MHC-II"
df2$Program[df2$group=="Cluster_13"] <- "Estrogen_Response"
df2$Program[df2$group=="Cluster_14"] <- "Epithelial_Senescence"
df2$Program[df2$group=="Cluster_15"] <- "Secreted"
df2$Program[df2$group=="Cluster_16"] <- "Chromatin_1"
df2$Program[df2$group=="Cluster_17"] <- "Chromatin_2"
df2$Program[df2$group=="Cluster_18"] <- "Hypoxia/Lipid_Metabolism"
df2$Program[df2$group=="Cluster_19"] <- "Translation_Initiation"
df2$Program[df2$group=="Cluster_20"] <- "Skin_Pigmentation"
df2$Program[df2$group=="Cluster_21"] <- "p53-Dependent_Senescence"
df2$Program <- factor(df2$Program, levels = c("CellCycle_G2M_1", "CellCycle_G2M_2", "CellCycle_G2M_3", 
                                              "CellCycle_G1S", "Proteasomal_Degradation", "Stress", "EMT-I", 
                                              "EMT-II_1", "EMT-II_2", "EMT-III", "EMT-IV", 
                                              "Interferon/MHC-II", "Estrogen_Response", "Epithelial_Senescence",
                                              "Secreted", "Chromatin_1", "Chromatin_2", "Hypoxia/Lipid_Metabolism", "Translation_Initiation", 
                                              "Skin_Pigmentation", "p53-Dependent_Senescence"))
p2 <- ggplot(df2, aes(x=x,y=y))+
  geom_tile(aes(fill=Program))+
  scale_x_continuous(expand = c(0,0))+
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "left",
        legend.title = element_blank(),
        legend.text = element_text(size = 15)) +
  scale_fill_manual(values = c(fig.colors, dark.fig.colors))
p1 <- ggplot(data = nmf_intersect_meltI_NEW, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +                                
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))
p3 <- p1 %>% insert_left(p2,width = 0.05)


# Style2
nmf_intersect_meltI_NEW <- reshape2::melt(nmf_intersect_original[rev(inds_sorted), rev(inds_sorted)])
nmf_intersect_meltI_NEW$group <- "Cluster_Other"
for (i in names(Cluster_list3)) {
  nmf_intersect_meltI_NEW$group[nmf_intersect_meltI_NEW$Var1 %in% Cluster_list3[[i]]] <- i
}
nmf_intersect_meltI_NEW$group <- factor(nmf_intersect_meltI_NEW$group, levels = mixedsort(unique(nmf_intersect_meltI_NEW$group)))
df2 <- tibble(y=nmf_intersect_meltI_NEW$Var1) %>% as.data.frame()
df2$x <- 1
df2$group <- nmf_intersect_meltI_NEW$group
df2$Program <- NA
unique(df2$group)
df2$Program[df2$group=="Cluster_1"] <- "Cell Cycle G2/M (1)"
df2$Program[df2$group=="Cluster_2"] <- "Cell Cycle G2/M (2)"
df2$Program[df2$group=="Cluster_3"] <- "Cell Cycle G2/M (3)"
df2$Program[df2$group=="Cluster_4"] <- "Cell Cycle G1/S"
df2$Program[df2$group=="Cluster_5"] <- "Proteasomal Degradation"
df2$Program[df2$group=="Cluster_6"] <- "Stress"
df2$Program[df2$group=="Cluster_7"] <- "EMT-I"
df2$Program[df2$group=="Cluster_8"] <- "EMT-II (1)"
df2$Program[df2$group=="Cluster_9"] <- "EMT-II (2)"
df2$Program[df2$group=="Cluster_10"] <- "EMT-III"
df2$Program[df2$group=="Cluster_11"] <- "EMT-IV"
df2$Program[df2$group=="Cluster_12"] <- "Interferon/MHC-II"
df2$Program[df2$group=="Cluster_13"] <- "Estrogen Response"
df2$Program[df2$group=="Cluster_14"] <- "Epithelial Senescence"
df2$Program[df2$group=="Cluster_15"] <- "Secreted"
df2$Program[df2$group=="Cluster_16"] <- "Chromatin (1)"
df2$Program[df2$group=="Cluster_17"] <- "Chromatin (2)"
df2$Program[df2$group=="Cluster_18"] <- "Hypoxia/Lipid Metabolism"
df2$Program[df2$group=="Cluster_19"] <- "Translation Initiation"
df2$Program[df2$group=="Cluster_20"] <- "Skin Pigmentation"
df2$Program[df2$group=="Cluster_21"] <- "p53-Dependent Senescence"
df2$Program <- factor(df2$Program, levels = c("Cell Cycle G2/M (1)", "Cell Cycle G2/M (2)", "Cell Cycle G2/M (3)",
                                              "Cell Cycle G1/S", "Proteasomal Degradation", "Stress", "EMT-I", 
                                              "EMT-II (1)","EMT-II (2)", "EMT-III", "EMT-IV", 
                                              "Interferon/MHC-II", "Estrogen Response", "Epithelial Senescence",
                                              "Secreted", "Chromatin (1)", "Chromatin (2)", "Hypoxia/Lipid Metabolism",
                                              "Translation Initiation", "Skin Pigmentation", "p53-Dependent Senescence"))

p2 <- ggplot(df2, aes(x=x,y=y))+
  geom_tile(aes(fill=Program))+
  scale_x_continuous(expand = c(0,0))+
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "left",
        legend.title = element_blank(),
        legend.text = element_text(size = 15)) +
  scale_fill_manual(values = c(fig.colors, dark.fig.colors))
p1 <- ggplot(data = nmf_intersect_meltI_NEW, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) +
  geom_tile() +
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))
p3 <- p1 %>% insert_left(p2,width = 0.05)


lapply(Cluster_list3, length) %>% unlist
total.len <- length(unlist(Cluster_list3)) + 0.5
total.len
len1 <- total.len-sum((lapply(Cluster_list3, length) %>% unlist)[1:3])
len2 <- len1-sum((lapply(Cluster_list3, length) %>% unlist)[4])
len3 <- len2-sum((lapply(Cluster_list3, length) %>% unlist)[5])
len4 <- len3-sum((lapply(Cluster_list3, length) %>% unlist)[6])
len5 <- len4-sum((lapply(Cluster_list3, length) %>% unlist)[7:11])
len6 <- len5-sum((lapply(Cluster_list3, length) %>% unlist)[12])
len7 <- len6-sum((lapply(Cluster_list3, length) %>% unlist)[13:14])
len8 <- len7-sum((lapply(Cluster_list3, length) %>% unlist)[15])
len9 <- len8-sum((lapply(Cluster_list3, length) %>% unlist)[16:17])
len10 <- len9-sum((lapply(Cluster_list3, length) %>% unlist)[18])
len11 <- len10-sum((lapply(Cluster_list3, length) %>% unlist)[19])
len12 <- len11-sum((lapply(Cluster_list3, length) %>% unlist)[20])
len13 <- len12-sum((lapply(Cluster_list3, length) %>% unlist)[21])

p2 <- ggplot(df2, aes(x=x,y=y))+
  geom_tile(aes(fill=Program))+
  scale_x_continuous(expand = c(0,0))+
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "left",
        legend.title = element_blank(),
        legend.text = element_text(size = 15)) +
  scale_fill_manual(values = c(fig.colors, dark.fig.colors))
p1 <- ggplot(data = nmf_intersect_meltI_NEW, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) +
  geom_tile() +
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1)) +
  geom_rect(aes(xmin = len1, xmax = total.len, ymin = len1, ymax = total.len), color = "black", fill = NA, linetype = "longdash", size = 0.7) + # Cluster123 307
  geom_rect(aes(xmin = len2, xmax = len1, ymin = len2, ymax = len1), color = "black", fill = NA, linetype = "longdash", size = 0.7) + # Cluster4
  geom_rect(aes(xmin = len3, xmax = len2, ymin = len3, ymax = len2), color = "black", fill = NA, linetype = "longdash", size = 0.7) + # Cluster5
  geom_rect(aes(xmin = len4, xmax = len3, ymin = len4, ymax = len3), color = "black", fill = NA, linetype = "longdash", size = 0.7) + # Cluster6
  geom_rect(aes(xmin = len5, xmax = len4, ymin = len5, ymax = len4), color = "black", fill = NA, linetype = "longdash", size = 0.7) + # Cluster7891011 185
  geom_rect(aes(xmin = len6, xmax = len5, ymin = len6, ymax = len5), color = "black", fill = NA, linetype = "longdash", size = 0.7) + # Cluster12
  geom_rect(aes(xmin = len7, xmax = len6, ymin = len7, ymax = len6), color = "black", fill = NA, linetype = "longdash", size = 0.7) + # Cluster1314 67
  geom_rect(aes(xmin = len8, xmax = len7, ymin = len8, ymax = len7), color = "black", fill = NA, linetype = "longdash", size = 0.7) + # Cluster15
  geom_rect(aes(xmin = len9, xmax = len8, ymin = len9, ymax = len8), color = "black", fill = NA, linetype = "longdash", size = 0.7) + # Cluster1617 35
  geom_rect(aes(xmin = len10, xmax = len9, ymin = len10, ymax = len9), color = "black", fill = NA, linetype = "longdash", size = 0.7) + # Cluster18
  geom_rect(aes(xmin = len11, xmax = len10, ymin = len11, ymax = len10), color = "black", fill = NA, linetype = "longdash", size = 0.7) + # Cluster19
  geom_rect(aes(xmin = len12, xmax = len11, ymin = len12, ymax = len11), color = "black", fill = NA, linetype = "longdash", size = 0.7) + # Cluster20
  geom_rect(aes(xmin = len13, xmax = len12, ymin = len13, ymax = len12), color = "black", fill = NA, linetype = "longdash", size = 0.7) # Cluster21
p3 <- p1 %>% insert_left(p2,width = 0.05)

# Style3
nmf_intersect_meltI_NEW <- reshape2::melt(nmf_intersect_original[rev(inds_sorted), rev(inds_sorted)])
nmf_intersect_meltI_NEW$group <- "Cluster_Other"
for (i in names(Cluster_list3)) {
  nmf_intersect_meltI_NEW$group[nmf_intersect_meltI_NEW$Var1 %in% Cluster_list3[[i]]] <- i
}
nmf_intersect_meltI_NEW$group <- factor(nmf_intersect_meltI_NEW$group, levels = mixedsort(unique(nmf_intersect_meltI_NEW$group)))
df2 <- tibble(y=nmf_intersect_meltI_NEW$Var1) %>% as.data.frame()
df2$x <- 1
df2$group <- nmf_intersect_meltI_NEW$group

df2$Program <- NA
unique(df2$group)
df2$Program[df2$group=="Cluster_1"] <- "CellCycle_G2M_1"
df2$Program[df2$group=="Cluster_2"] <- "CellCycle_G2M_2"
df2$Program[df2$group=="Cluster_3"] <- "CellCycle_G2M_3"
df2$Program[df2$group=="Cluster_4"] <- "CellCycle_G1S"
df2$Program[df2$group=="Cluster_5"] <- "Proteasomal_Degradation"
df2$Program[df2$group=="Cluster_6"] <- "Stress"
df2$Program[df2$group=="Cluster_7"] <- "EMT_I"
df2$Program[df2$group=="Cluster_8"] <- "EMT_II_1"
df2$Program[df2$group=="Cluster_9"] <- "EMT_II_2"
df2$Program[df2$group=="Cluster_10"] <- "EMT_III"
df2$Program[df2$group=="Cluster_11"] <- "EMT_IV"
df2$Program[df2$group=="Cluster_12"] <- "Interferon__MHC_II"
df2$Program[df2$group=="Cluster_13"] <- "Estrogen_Response"
df2$Program[df2$group=="Cluster_14"] <- "Epithelial_Senescence"
df2$Program[df2$group=="Cluster_15"] <- "Secreted"
df2$Program[df2$group=="Cluster_16"] <- "Chromatin_1"
df2$Program[df2$group=="Cluster_17"] <- "Chromatin_2"
df2$Program[df2$group=="Cluster_18"] <- "Hypoxia__Lipid_Metabolism"
df2$Program[df2$group=="Cluster_19"] <- "Translation_Initiation"
df2$Program[df2$group=="Cluster_20"] <- "Skin_Pigmentation"
df2$Program[df2$group=="Cluster_21"] <- "p53_Dependent_Senescence"
df2$Program <- factor(df2$Program, levels = c("CellCycle_G2M_1", "CellCycle_G2M_2", "CellCycle_G2M_3", 
                                              "CellCycle_G1S", "Proteasomal_Degradation", "Stress", "EMT_I", 
                                              "EMT_II_1", "EMT_II_2", "EMT_III", "EMT_IV", 
                                              "Interferon__MHC_II", "Estrogen_Response", "Epithelial_Senescence",
                                              "Secreted", "Chromatin_1", "Chromatin_2", "Hypoxia__Lipid_Metabolism", "Translation_Initiation", 
                                              "Skin_Pigmentation", "p53_Dependent_Senescence"))
p2 <- ggplot(df2, aes(x=x,y=y))+
  geom_tile(aes(fill=Program))+
  scale_x_continuous(expand = c(0,0))+
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "left",
        legend.title = element_blank(),
        legend.text = element_text(size = 15)) +
  scale_fill_manual(values = c(fig.colors, dark.fig.colors))
p1 <- ggplot(data = nmf_intersect_meltI_NEW, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +                                
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))
p3 <- p1 %>% insert_left(p2,width = 0.05)




MP_list4 <- MP_list3
names(MP_list4) <- levels(df2$Program)
MP_list4 %>% unlist %>% as.character
MP_list4.df <- as.data.frame(MP_list4)
MP_list4.df <- read.csv("E:/Master_2nd_year/CellLine/4.Main_results/CellLine_Atlas_outs3.0/Part2/8.InputFeatureForCapsNet_MP_list4.csv")
names(MP_list4.df)
MP_list4.df$Translation.Initiation

intersect(MP_list4.df$Translation.Initiation, Nat.MP.df$gene[Nat.MP.df$Nat.MP=="MP11 Translation initiation"])



#### 5.13 tissue ####

## CCLE
test <- readxl::read_xlsx("./IF38.33.NG.41588_2020_726_MOESM3_ESM.xlsx", skip=3)
test$CellLine_name <- lapply(test$Name, function(x){strsplit(x, "_")[[1]][1]}) %>% unlist()
test

organ.df <- tibble(CellLine_name="", Cancer_type="")
for (i in 1:nrow(test)) {
  test.df <- tibble(CellLine_name=test$CellLine_name[i], Cancer_type=test$Primary_Disease[i])
  organ.df <- rbind(organ.df, test.df)
}
organ.df <- organ.df[-1,]
table(organ.df$Cancer_type)

## CL
CL1.name <- c("BT474", "HCC1937", "MCF10A", "MCF7", "MDAMB231", "SKBR3", "T47D",
              "A549", "MRC5", "NCIH446", "NCIH460", "PC9")
CL2.name <- c("HELA","SIHA",
              "786O","ACHN",
              "HUH7","PLCPRF5","SKHEP1",
              "EC109","KYSE150","KYSE510","TE1",
              "ASPC1","CAPAN2","MIAPACA2","PANC1","SW1990",
              "BCPAP",
              "T24")
CL3.name <- c("C33A",
              "AN3CA","HEC1A","HEC1B",
              "HCT116","HT29","RKO","SW620",
              "KYSE410",
              "A2780","SKOV3",
              "DU145","LNCaP","PC3",
              "A431",
              "CASKI",
              "DLD1",
              "AGS","NCIN87",
              "J82")

CL1.name.cancer <- c("Breast Cancer", "Breast Cancer", "Breast Epithelium", "Breast Cancer", "Breast Cancer", "Breast Cancer", "Breast Cancer", 
                     "Lung Cancer", "Fibroblast", "Lung Cancer", "Lung Cancer", "Lung Cancer")
test.df <- tibble(CellLine_name=CL1.name, Cancer_type=CL1.name.cancer)
organ.df <- rbind(organ.df, test.df)
CL2.name.cancer <- c("Cervical Cancer", "Cervical Cancer", 
                     "Kidney Cancer", "Kidney Cancer", 
                     "Liver Cancer", "Liver Cancer", "Liver Cancer", 
                     "Esophageal Cancer", "Esophageal Cancer", "Esophageal Cancer", "Esophageal Cancer", 
                     "Pancreatic Cancer", "Pancreatic Cancer", "Pancreatic Cancer", "Pancreatic Cancer", "Pancreatic Cancer",
                     "Thyroid Cancer",
                     "Bladder Cancer")
test.df <- tibble(CellLine_name=CL2.name, Cancer_type=CL2.name.cancer)
organ.df <- rbind(organ.df, test.df)
CL3.name.cancer <- c("Cervical Cancer", 
                     "Endometrial/Uterine Cancer", "Endometrial/Uterine Cancer", "Endometrial/Uterine Cancer", 
                     "Colon/Colorectal Cancer", "Colon/Colorectal Cancer", "Colon/Colorectal Cancer", "Colon/Colorectal Cancer", 
                     "Esophageal Cancer", 
                     "Ovarian Cancer", "Ovarian Cancer", 
                     "Prostate Cancer", "Prostate Cancer", "Prostate Cancer", 
                     "Skin Cancer", 
                     "Cervical Cancer", 
                     "Colon/Colorectal Cancer", 
                     "Gastric Cancer", "Gastric Cancer",
                     "Bladder Cancer")
test.df <- tibble(CellLine_name=CL3.name, Cancer_type=CL3.name.cancer)
organ.df <- rbind(organ.df, test.df)


## BC7
mtdt <- read.csv(paste0(public_dir, "GSM4285803/GSM4285803_scRNA_metaInfo.csv"))
mtdt <- column_to_rownames(mtdt, "X")
head(mtdt);dim(mtdt)

BC7.name <- unique(mtdt$CellType)
BC7.name <- BC7.name[-which(BC7.name %in% organ.df$CellLine_name)]
BC7.name <- BC7.name[-c(2,3)]
BC7.name
BC7.name.cancer <- rep("Breast Cancer", length(BC7.name))
test.df <- tibble(CellLine_name=BC7.name, Cancer_type=BC7.name.cancer)
organ.df <- rbind(organ.df, test.df)


## BC32
BC32 <- readRDS(paste0(public_dir, "GSE173634/RAW.UMI.counts.BC.cell.lines.rds"))
ids=annoGene(rownames(BC32),'ENSEMBL','human')
head(ids)
ids=ids[!duplicated(ids$SYMBOL),]
ids=ids[!duplicated(ids$ENSEMBL),]
dim(ids)
kp = rownames(BC32) %in% ids$ENSEMBL
table(kp)
BC32 <- BC32[kp,]
rownames(BC32)=ids$SYMBOL[match( rownames(BC32) , ids$ENSEMBL)]
BC32[1:4,1:4]
BC32[1:4,(ncol(BC32)-4):ncol(BC32)]

BC32.seu <- CreateSeuratObject(counts = BC32, min.cells = 3, min.features = 200)
BC32.seu$CellLine <- lapply(colnames(BC32.seu), function(x){
  strsplit(x, "_")[[1]][1]
}) %>% unlist()
table(BC32.seu$CellLine)

BC32.name <- unique(BC32.seu$CellLine)
BC32.name <- BC32.name[-which(BC32.name %in% organ.df$CellLine_name)]
BC32.name
BC32.name.cancer <- rep("Breast Cancer", length(BC32.name))
test.df <- tibble(CellLine_name=BC32.name, Cancer_type=BC32.name.cancer)
organ.df <- rbind(organ.df, test.df)


## LCL
LCL.name <- c("LCL777B958", "LCL461B958", "LCL777B958M81AGGR")
LCL.name.cancer <- rep("B lymphoblastoma", length(LCL.name))
test.df <- tibble(CellLine_name=LCL.name, Cancer_type=LCL.name.cancer)
organ.df <- rbind(organ.df, test.df)
organ.df

organ.df <- read.csv("./7.organ.df.csv")





#### 5.14 phyper & IS score ####

cal_CancerType <- function(organ.df=organ.df, 
                           Cluster_list=Cluster_list3,
                           Cluster_idx=15) {
  res.df <- organ.df$Cancer_type[organ.df$CellLine_name %in%
                                   c(lapply(Cluster_list[[paste0("Cluster_", Cluster_idx)]], function(x){strsplit(x, "_")[[1]][2]}) %>%
                                       unlist %>% unique)] %>% table %>% sort(decreasing = T)

  phyper.df <- tibble(Program.name="", pvalue=0)
  for (i in 1:length(res.df)) {
    N = nrow(organ.df)
    M = table(organ.df$Cancer_type)[names(res.df[i])][[1]]
    n = sum(res.df)
    m = res.df[i][[1]]
    
    pvalue = phyper(m,M,(N-M),n,lower.tail=F)
    phyper.df <- rbind(phyper.df, tibble(Program.name=names(res.df[i]), 
                                         pvalue=pvalue))
  }
  phyper.df <- phyper.df[-1,]
  if (nrow(phyper.df)!=1) {
    phyper.df$padjust <- p.adjust(phyper.df$pvalue, method='fdr')
  } else {
    phyper.df$padjust <- phyper.df$pvalue
  }
  
  res.df <- as.data.frame(res.df)
  if (length(res.df)==1) {
    res.df$ClusterFreq <- res.df$res.df
    res.df$res.df <- rownames(res.df)
    colnames(res.df)[1] <- "CancerType"
  } else {
    colnames(res.df)[1] <- "CancerType"
    colnames(res.df)[2] <- "ClusterFreq"
  }
  test2 <- as.data.frame(table(organ.df$Cancer_type))
  colnames(test2)[1] <- "CancerType"
  res.df <- left_join(res.df, y = test2)
  colnames(res.df)[3] <- "OrganFreq"
  
  res.df$inter.score <- res.df$ClusterFreq / res.df$OrganFreq
  res.df$inter.score <- res.df$inter.score / sum(res.df$inter.score)
  
  return(list(phyper.df, res.df))
}

organ.phyper.df <- tibble(MP="", CancerType="", ClusterFreq=0, OrganFreq=0, inter.score=0, pvalue=0, padjust=0)
for (i in 1:length(Cluster_list3)) {
  test <- cal_CancerType(organ.df=organ.df, 
                         Cluster_list=Cluster_list3, 
                         Cluster_idx=i)
  colnames(test[[1]])[1] <- "CancerType"
  test <- right_join(test[[1]] %>% arrange(padjust), test[[2]])
  test$MP <- names(Cluster_list3)[i]
  test <- test[,c(7, 1, 4, 5, 6, 2, 3)]
  organ.phyper.df <- rbind(organ.phyper.df, test)
}
organ.phyper.df <- organ.phyper.df[-1,]
organ.phyper.df


df3 <- df2
head(df3)
colnames(df3)[3] <- 'MP'
df3 <- df3[,c(3,4)]
df3 <- df3[!duplicated.data.frame(df3),]
df3

organ.phyper.df <- left_join(x = organ.phyper.df, y = df3, by = 'MP')
colnames(organ.phyper.df)
head(organ.phyper.df)


organ.phyper.df2$padjust[organ.phyper.df2$padjust==0] <- min(organ.phyper.df2$padjust[organ.phyper.df2$padjust!=0])

ggplot(organ.phyper.df2, aes(x=Program, y=CancerType, color=-log(padjust), size=inter.score)) +
  geom_point() +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(color = 'black', size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  ) +
  scale_color_continuous(low = fig.colors[11], high = fig.colors[1]) +
  labs(size='Score')









