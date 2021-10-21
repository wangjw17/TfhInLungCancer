####### Analysis of data from Kim et al. (Nat Comm 2020) ######
### Jiawei Wang, Zhao lab, 2020/10
### Here below is the code used to analyze single-cell RNA-seq data of T cells and B cells from
### Kim et al. R package Seurat was used for data input, cell type selection, cell clustering,
### marker gene expression calculation and visualization.

###### Clustering of T/NK cells from LUNG + LN ######
### Extraction of LUNG + LN cells ####

library(ggplot2)
library(dplyr)
library(Seurat)
library(stringr)
library(data.table)

system('head -1 GSE131907_Lung_Cancer_raw_UMI_matrix.txt > header_count.txt')
umi_header <- readLines('header_count.txt')
samples <- strsplit(umi_header, split="\t")[[1]]
meta_umi <- data.frame(raw=samples[-1],
                       barcode = samples[-1] %>% strsplit(split="_") %>% sapply("[",1),
                       source = samples[-1] %>% strsplit(split="_") %>% sapply("[",2),
                       index = samples[-1] %>% strsplit(split="_") %>% sapply("[",3))

mat <- fread('GSE131907_Lung_Cancer_raw_UMI_matrix.txt')
dim(mat)

annot <- read.table('GSE131907_Lung_Cancer_cell_annotation.txt', sep="\t", header=T)
dim(annot)
head(annot)
sum(annot$Sample_Origin %in% c("mLN","nLN","nLung","tLung"))
annot$number <- annot$Index %>% as.character %>% strsplit(split="_") %>% sapply("[",3)
annot <- annot[match(meta_umi$raw, annot$Index),]

# sum(annot$Index == meta_umi$raw)
# sum(meta_umi$raw == names(mat)[-1])

sel0 <- intersect(which(annot$Sample_Origin %in% c("nLung","tLung","mLN","nLN")), which(annot$Cell_type.refined == "T/NK cells"))
annot_ll <- annot[sel0,]
# sum(annot$Sample_Origin  %in% c("nLung","tLung","mLN","nLN") & annot$Cell_type.refined == "T/NK cells", na.rm=T)

class(mat) <- class(as.data.frame(mat))
mat_ll <- mat[,1+sel0]
# sum(names(mat_ll) == annot_ll$Index)

genes <- mat$Index

annot_ll$index_new <- NA
for (num in levels(factor(annot_ll$number))){
    id <- which(annot_ll$number==num)
    annot_ll$index_new[id] <- paste0(num, "_", 1:length(id))
}

names(mat_ll) <- annot_ll$index_new
save(annot_ll, genes, mat_ll, file="kim_umi_ll.RData")





### Clustering of LUNG/LN T/NK cells with Seurat ####

rm(list=ls())
load('kim_umi_ll.RData')
# rm(annot, cel0, id, lungcell, mat, meta_umi, num, samples, sel0, umi_header)

rownames(mat_ll) <- genes
kp <- CreateSeuratObject(counts = mat_ll, project = "KimDataset", min.cells = 3, min.features = 200)
summary(kp)

kp[["percent.mt"]] <- PercentageFeatureSet(kp, pattern = "^MT-")

VlnPlot(kp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)

plot1 <- FeatureScatter(kp, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(kp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# CombinePlots(plots = list(plot1, plot2))
# kp <- subset(kp, subset = nCount_RNA < 15000 & nFeature_RNA < 5000 & percent.mt < 10)

## normalization
kp <- NormalizeData(kp, normalization.method = "LogNormalize", scale.factor = 10000)
# kp <- NormalizeData(kp, normalization.method = "RC", scale.factor = 1e6) ## for heatmap of relative expression
kp <- FindVariableFeatures(kp, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(kp), 10)
plot1 <- VariableFeaturePlot(kp)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
# CombinePlots(plots = list(plot1, plot2))

## dimension reduction
all.genes <- rownames(kp)
kp <- ScaleData(kp, features = all.genes)
kp <- RunPCA(kp, features = VariableFeatures(kp))
print(kp[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(kp, dims = 1:5, reduction = "pca")
DimHeatmap(kp, dims = 1:6, cells = 500, balanced = TRUE)
DimPlot(kp, reduction="pca")

## clustering cells
kp <- FindNeighbors(kp, dims = 1:20)
kp <- FindClusters(kp, resolution = 0.5)
# head(Idents(kp), 5)
table(Idents(kp))

## UMPA
# kp <- RunUMAP(kp, dims = 1:15)
DimPlot(kp, reduction = "umap")
# kp <- RunTSNE(kp, dims = 1:20)
# DimPlot(kp, reduction = "tsne")

rm(list=ls())
load('kim_umi_ll_seurat.RData')

pdf('umap_kim_ll.pdf', width=7, height=6)
DimPlot(kp, reduction = "umap")#, plot.title="Clutering of nLung/tLung cells, original clusters")
dev.off()

save(kp, file="kim_umi_ll_seurat.RData")



### Marker gene inspection #### 

load('kim_umi_ll_seurat.RData')
load('kim_umi_ll.RData')
rm(mat_ll)
annot_ll$cluster <- Idents(kp)

markers.t <- c("CD3D", "CD3E", "CD3G")
markers.tfh <- "IL21, IL4, CXCL13, CXCR5, TOX, TOX2, BCL6, PDCD1, ICOS, CD200, MAF, ASCL2, SH2D1A, SLAMF6, BTLA, GNG4, POU2AF1, IL6ST" %>% strsplit(split=", ") %>% unlist
markers.th17 <- "IL17A, IL17F, RORA, RORC, STAT3, NR4A2" %>% strsplit(split=", ") %>% unlist
markers.th1 <- "IFNG, IL2, TBX21, STAT4, IL12RB2, BHLHE40, CXCR3" %>% strsplit(split=", ") %>% unlist
markers.th2 <- "IL4, IL5, IL13, GATA3" %>% strsplit(split=", ") %>% unlist
markers.treg <- "IL2RA, FOXP3, IKZF2, TNFRSF4, TNFRSF18" %>% strsplit(split=", ") %>% unlist

markers.nk <- "KLRC1, KLRD1, KLRF1" %>% strsplit(split=", ") %>% unlist #NKG2A=KLRC1
markers.naive <- "LEF1, TCF7, SELL, CCR7, CD8A, CD8B" %>% strsplit(split=", ") %>% unlist #TCF1=TCF7
markers.ctoxic <- "GZMA, GZMB, GZMK, PRF1, GNLY, IFNG, IL2, KLRG1, ZEB2, CX3CR1" %>% strsplit(split=", ") %>% unlist
markers.inhr <- "LAG3, TIGIT, PDCD1, HAVCR2, CTLA4" %>% strsplit(split=", ") %>% unlist


pdf('tsne_kim_ll.pdf')
DimPlot(kp, reduction = "tsne")#, plot.title="Clutering of nLung/tLung cells, original clusters")
dev.off()

pdf('marker_kim_ll_tnk.pdf', width=20, height=12)
markers.tnk <- c(markers.t, markers.nk)
DimPlot(kp, reduction = "umap")
VlnPlot(kp, features = markers.tnk, pt.size=0)
FeaturePlot(kp, features = markers.tnk)
dev.off()

pdf('marker_kim_ll_tcd4.pdf', width=20, height=12)
VlnPlot(kp, features = markers.tfh, pt.size=0)
VlnPlot(kp, features = markers.th17, pt.size=0)
VlnPlot(kp, features = markers.th1, pt.size=0)
VlnPlot(kp, features = markers.th2, pt.size=0)
VlnPlot(kp, features = markers.treg, pt.size=0)
dev.off()

pdf('marker_kim_ll_nk.pdf', width=20, height=12)
VlnPlot(kp, features = markers.naive, pt.size=0)
VlnPlot(kp, features = markers.ctoxic, pt.size=0)
VlnPlot(kp, features = markers.inhr, pt.size=0)
dev.off()

##calculate differential expression
kp.markers <- FindAllMarkers(kp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top30 <- kp.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)

top30[top30$cluster==11,]

write.csv(top30, file="top30.csv")

## COMMENT
# CD8: 2,4,6,9,12
# 
# NK: 3,13
# 
# CD4: 0(naive),1(some Th),5(Treg),7(?),8(naive),10(Tfh),11(Th17),14(Treg)

### Plot heatmap #### 

library(ggplot2)
library(dplyr)
library(Seurat)
library(stringr)
library(data.table)

rm(list=ls())
load('kim_umi_ll_seurat.RData')
clusters <- Idents(kp)

## Instead, construct Seurat object with NormalizeData(method="RC")
load('kim_umi_ll.RData')
rownames(mat_ll) <- genes
kp <- CreateSeuratObject(counts = mat_ll, project = "KimDataset", min.cells = 3, min.features = 200)
kp <- NormalizeData(kp, normalization.method = "RC", scale.factor = 1e6) ## for heatmap of relative expression
Idents(kp) <- clusters

save(kp, file="kim_umi_ll_seurat_rc.RData")

markers.tfh <- "IL21, IL4, CXCL13, CXCR5, TOX, TOX2, BCL6, PDCD1, ICOS, CD200, MAF, ASCL2, SH2D1A, SLAMF6, BTLA" %>% strsplit(split=", ") %>% unlist
markers.th17 <- "IL17A, IL17F, RORA, RORC, STAT3, NR4A2" %>% strsplit(split=", ") %>% unlist
markers.th1 <- "IFNG, IL2, TBX21, STAT4, IL12RB2, BHLHE40, CXCR3" %>% strsplit(split=", ") %>% unlist
markers.th2 <- "IL4, IL5, IL13, GATA3" %>% strsplit(split=", ") %>% unlist
markers.treg <- "IL2RA, FOXP3, IKZF2, TNFRSF4, TNFRSF18" %>% strsplit(split=", ") %>% unlist
markers.t <- "CD3D, CD3E, CD3G, CD8A, CD8B, CD4" %>% strsplit(split=", ") %>% unlist
markers.nk <- "KLRC1, KLRD1, KLRF1" %>% strsplit(split=", ") %>% unlist #NKG2A=KLRC1
markers.naive <- "LEF1, TCF7, SELL, CCR7" %>% strsplit(split=", ") %>% unlist #TCF1=TCF7
markers.ctoxic <- "GZMA, GZMB, GZMK, PRF1, GNLY, IFNG, IL2, KLRG1, ZEB2, CX3CR1" %>% strsplit(split=", ") %>% unlist
markers.inhr <- "LAG3, TIGIT, PDCD1, HAVCR2, CTLA4" %>% strsplit(split=", ") %>% unlist
markers.all <- c(markers.tfh, markers.th17, markers.th1, markers.th2, markers.treg, markers.t, markers.nk, markers.naive, markers.ctoxic, markers.inhr) %>% unique
table(markers.all) ##IFNG, PDCD1, IL4, IL2

markers.tfh <- "IL21, IL4, CXCR5, CXCL13, BCL6, TOX, TOX2, BTLA, PDCD1, CD200, MAF, SH2D1A, ASCL2" %>% strsplit(split=", ") %>% unlist
markers.treg <- "IL2RA, FOXP3, IKZF2" %>% strsplit(split=", ") %>% unlist
markers.th17 <- "IL17A, IL17F, RORA, RORC" %>% strsplit(split=", ") %>% unlist
markers.th1 <- "IFNG, TBX21, CXCR3, IL12RB2" %>% strsplit(split=", ") %>% unlist
markers.th2 <- "IL5, IL13, GATA3" %>% strsplit(split=", ") %>% unlist
markers.naive <- "LEF1, TCF7, SELL, CCR7" %>% strsplit(split=", ") %>% unlist #TCF1=TCF7
markers.sel <-  c(markers.tfh, markers.treg, markers.th17, markers.th1, markers.th2, markers.naive)

kp.mk <- kp[['RNA']][markers.sel,] %>% as.matrix %>% t()
kp.mk <- log(base=2, kp.mk+1) %>% as.data.frame
kp.mk$cluster <- Idents(kp)
kp.ag <- aggregate(kp.mk[,1:(ncol(kp.mk)-1)], list(kp.mk$cluster), mean)

# kp.ag <- apply(kp.ag[,-1], 1, function(x) scale(x, scale = F, center=T))
kp.ag <- scale(kp.ag[,-1])
rownames(kp.ag) <- levels(kp.mk$cluster)

## Plotting order of clusters: 14(pTfh), 10(Tfh), 5(Treg), 1(Th17,Th1,Th2), 7(Th17,naive), 0, 8(naive)

clusters.sel <- c(14,10,5,1,7,0,8,11)
df.ag <- melt(kp.ag)
df.ag <- df.ag[df.ag$Var1 %in% clusters.sel,]
names(df.ag) <- c("Cluster","Gene","RelAb")
range(df.ag$RelAb)
df.ag$Gene <- factor(df.ag$Gene, levels=markers.sel %>% rev)
df.ag$Cluster <- factor(df.ag$Cluster, levels=clusters.sel)

p1 <- ggplot(df.ag, aes(x = Cluster, y=Gene, fill=RelAb)) +
  geom_tile(colour="white", size=.5) +
  scale_fill_gradient2(low = "blue4", high = "red3", mid = "white",# midpoint = 0, limit = c(-4,4), 
                       space = "Lab", name="Relative\nAbundance") +
  scale_y_discrete(position = "right") +
  labs(x=NULL,y=NULL)+
  theme(axis.text.y.right = element_text(face="italic"),
        panel.background=element_blank())
# pdf('heatmap_kim_ll_cd4.pdf', height=4, width=4) ##(3): 4x4
p1
# dev.off()

table(annot_ll$Cell_subtype, Idents(kp)) %>% write.csv(file="celltype_cluster_assignment.csv")

### Experiment with ScaleData, LogNormalize and relative abundance measures ####
## The idea of relative abundance comes from https://www.nature.com/articles/nature24489 but 
## depicted differently in Kim et al.

#use CD200 as an example
names(kp.mk)[10]
hist(kp.mk[,10])

load('kim_umi_ll.RData')

## repeat with previous steps
cd200.raw <- kp[['RNA']]['CD200',] %>% as.numeric
## then logNormalize
cd200.log <- kp[['RNA']]['CD200',] %>% as.numeric
## then scale
cd200.scl <- kp[['RNA']]['CD200',] %>% as.numeric

hist(kp@assays$RNA['CD200',] %>% as.numeric)

kp <- ScaleData(kp, do.scale = F, do.center = T,assay='RNA', features=rownames(kp))

cd200 <- kp[['RNA']]['CD200',] %>% as.numeric
hist(cd200[cd200>0])

## COMMENT 
# Conclusion: ScaleData doesn't perform on RNA assay. Still not sure why ScaleData doesn't take 
# effect on the data. But will continue munually with the procedures.

### Stack barplot of clusters ####
origin <- table(annot_ll$Sample_Origin, Idents(kp)) %>% as.data.frame
names(origin) <- c("Origin","Cluster","Number")
origin <- origin[origin$Cluster %in% clusters.sel & origin$Origin %in% c("mLN",'nLN','tLung','nLung'),]
origin$Cluster <- factor(origin$Cluster, levels=clusters.sel)
origin$Origin <- factor(origin$Origin, levels=c("nLung",'tLung','nLN','mLN'))

for(c in levels(origin$Cluster)) origin$Number[origin$Cluster==c] <- origin$Number[origin$Cluster==c]/sum(origin$Number[origin$Cluster==c])*100

p2 <- ggplot(origin, aes(x=Cluster, y=Number, fill=Origin)) +
    geom_bar(stat="identity", width=0.8) + 
    scale_fill_manual(values=c("gold1","firebrick","darkseagreen2","green4")) +
    scale_y_continuous(position = "right") +
    labs(x=NULL, y="Percentage", size=.5)+
    theme(panel.background=element_blank())

pdf('barplot_kim_ll_cd4_stack.pdf', width=6, height=2.5)
p2
dev.off()

library(grid)
library(gridExtra)

p1 <- ggplot(df.ag, aes(x = Cluster, y=Gene, fill=RelAb)) +
  geom_tile(colour="white", size=.5) +
  scale_fill_gradient2(low = "blue4", high = "red3", mid = "white", midpoint = 0, limit = c(-2,4), 
                       space = "Lab", name="Relative\nAbundance") +
  scale_y_discrete(position = "right") +
  labs(x=NULL,y=NULL)+
  theme(axis.text.y.right = element_text(face="italic"),
        panel.background=element_blank())

pdf('heatmap_kim_ll_cd4_stack.pdf', width=6, height=12)
grid.arrange(p2, p1, nrow = 2, widths = 30, heights = c(10, 80))
dev.off()






###### Analysis of B cells ######

### Cell selection ####
library(ggplot2)
library(dplyr)
library(Seurat)
library(stringr)
library(data.table)

rm(list=ls())
system('head -1 GSE131907_Lung_Cancer_raw_UMI_matrix.txt > header_count.txt')
umi_header <- readLines('header_count.txt')
samples <- strsplit(umi_header, split="\t")[[1]]
meta_umi <- data.frame(raw=samples[-1],
                       barcode = samples[-1] %>% strsplit(split="_") %>% sapply("[",1),
                       source = samples[-1] %>% strsplit(split="_") %>% sapply("[",2),
                       index = samples[-1] %>% strsplit(split="_") %>% sapply("[",3))

mat <- fread('GSE131907_Lung_Cancer_raw_UMI_matrix.txt')
dim(mat)

annot <- read.table('GSE131907_Lung_Cancer_cell_annotation.txt', sep="\t", header=T)
dim(annot)
head(annot)
sum(annot$Sample_Origin %in% c("mLN","nLN","nLung","tLung"))
annot$number <- annot$Index %>% as.character %>% strsplit(split="_") %>% sapply("[",3)
annot <- annot[match(meta_umi$raw, annot$Index),]

sel0 <- intersect(which(annot$Sample_Origin %in% c("nLung","tLung","mLN","nLN")), which(annot$Cell_type.refined == "B lymphocytes"))
annot_ll <- annot[sel0,]
# sum(annot$Sample_Origin  %in% c("nLung","tLung","mLN","nLN") & annot$Cell_type.refined == "T/NK cells", na.rm=T)

class(mat) <- class(as.data.frame(mat))
mat_ll <- mat[,1+sel0]
# sum(names(mat_ll) == annot_ll$Index)

genes <- mat$Index

annot_ll$index_new <- NA
for (num in levels(factor(annot_ll$number))){
    id <- which(annot_ll$number==num)
    annot_ll$index_new[id] <- paste0(num, "_", 1:length(id))
}

names(mat_ll) <- annot_ll$index_new
save(annot_ll, genes, mat_ll, file="kim_umi_llb.RData")



### Clustering of B cells ####

rownames(mat_ll) <- genes
kp <- CreateSeuratObject(counts = mat_ll, project = "KimDataset", min.cells = 3, min.features = 200)
summary(kp)

kp[["percent.mt"]] <- PercentageFeatureSet(kp, pattern = "^MT-")
VlnPlot(kp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)

plot1 <- FeatureScatter(kp, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(kp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# CombinePlots(plots = list(plot1, plot2))
# kp <- subset(kp, subset = nCount_RNA < 15000 & nFeature_RNA < 5000 & percent.mt < 10)

## normalization
kp <- NormalizeData(kp, normalization.method = "LogNormalize", scale.factor = 10000)
# kp <- NormalizeData(kp, normalization.method = "RC", scale.factor = 1e6) ## for heatmap of relative expression
kp <- FindVariableFeatures(kp, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(kp), 10)
plot1 <- VariableFeaturePlot(kp)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
# CombinePlots(plots = list(plot1, plot2))

## dimension reduction
all.genes <- rownames(kp)
kp <- ScaleData(kp, features = all.genes)
kp <- RunPCA(kp, features = VariableFeatures(kp))
print(kp[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(kp, dims = 1:5, reduction = "pca")
DimHeatmap(kp, dims = 1:6, cells = 500, balanced = TRUE)
DimPlot(kp, reduction="pca")

## clustering cells
kp <- FindNeighbors(kp, dims = 1:20)
kp <- FindClusters(kp, resolution = 0.25)
# head(Idents(kp), 5)
table(Idents(kp))

# clusters50 <- Idents(kp)
clusters25 <- Idents(kp)

## UMPA
kp <- RunUMAP(kp, dims = 1:15)
DimPlot(kp, reduction = "umap")
# kp <- RunTSNE(kp, dims = 1:20)
# DimPlot(kp, reduction = "tsne")

load('kim_umi_llb_seurat.RData')

pdf('umap_kim_llb_25.pdf', width=7, height=6)
DimPlot(kp, reduction = "umap")
dev.off()

# save(kp, annot_ll, clusters25, clusters50, file="kim_umi_llb_seurat.RData")
# load('kim_umi_llb_seurat.RData')

##calculate differential expression
kp.markers <- FindAllMarkers(kp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top30 <- kp.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)

write.csv(top30, file="top30_llb.csv")

assignment <- table(annot_ll$Cell_subtype, Idents(kp))
write.csv(assignment, file="celltype_cluster_assignment_b.csv")



### Heatmap of B cell markers (2 resolutions 0.25/0.5 ####

rm(list=ls())
# rownames(mat_ll) <- genes
# kp <- CreateSeuratObject(counts = mat_ll, project = "KimDataset", min.cells = 3, min.features = 200)
# kp <- NormalizeData(kp, normalization.method = "RC", scale.factor = 1e6) ## for heatmap of relative expression
# save(kp, clusters25, clusters50, annot_ll, file="kim_umi_llb_seurat_rc.RData")
load('kim_umi_llb_seurat_rc.RData')

markers.b1 <- "AICDA, BCL6, POLH, P2RY8, SEMA4A, FOXO1, BACH2, BATF, CD86, IRF8, DOCK8" %>% strsplit(split=", ") %>% unlist
markers.b2 <- "PRDM1, SDC1, XBP1, MZB1, JCHAIN, IGHA1, IGHA2, IGHG1, IGHG2" %>% strsplit(split=", ") %>% unlist
markers.b3 <- "MS4A1, HLA-DRA" %>% strsplit(split=", ") %>% unlist
markers.b4 <- "CXCR5, PAX5" %>% strsplit(split=", ") %>% unlist #"GZMB, IRF4" %>% strsplit(split=", ") %>% unlist

markers.sel <- c(markers.b1, markers.b2, markers.b3, markers.b4) %>% unique
Idents(kp) <- clusters25

kp.mk <- kp[['RNA']][markers.sel,] %>% as.matrix %>% t()
kp.mk <- log(base=2, kp.mk+1) %>% as.data.frame
kp.mk$cluster <- Idents(kp)
kp.ag <- aggregate(kp.mk[,1:(ncol(kp.mk)-1)], list(kp.mk$cluster), mean)

# kp.ag <- apply(kp.ag[,-1], 1, function(x) scale(x, scale = F, center=T))
kp.ag <- scale(kp.ag[,-1])
rownames(kp.ag) <- levels(kp.mk$cluster)
clusters.sel <- c(5, 4, 6, 0, 1, 2, 3, 8, 7)
df.ag <- melt(kp.ag)
df.ag <- df.ag[df.ag$Var1 %in% clusters.sel,]
names(df.ag) <- c("Cluster","Gene","RelAb")
range(df.ag$RelAb)
df.ag$Gene <- factor(df.ag$Gene, levels=markers.sel %>% rev)
df.ag$Cluster <- factor(df.ag$Cluster, levels=clusters.sel)

p1 <- ggplot(df.ag, aes(x = Cluster, y=Gene, fill=RelAb)) +
  geom_tile(colour="white", size=.5) +
  scale_fill_gradient2(low = "blue4", high = "red3", mid = "white", midpoint = 0, limit = c(-3,3), 
                       space = "Lab", name="Relative\nAbundance") +
  scale_y_discrete(position = "right") +
  labs(x=NULL,y=NULL)+
  theme(axis.text.y.right = element_text(face="italic"),
        panel.background=element_blank())
pdf('heatmap_kim_llb_b25.pdf', height=4, width=5.1) ##(3): 4x4
p1
dev.off()

### Stacked barplot for B cells ####
clusters <- clusters25
clusters.sel <- clusters.sel
origin <- table(annot_ll$Sample_Origin, clusters) %>% as.data.frame
names(origin) <- c("Origin","Cluster","Number")
origin <- origin[origin$Cluster %in% clusters.sel & origin$Origin %in% c("mLN",'nLN','tLung','nLung'),]
origin$Cluster <- factor(origin$Cluster, levels=clusters.sel)
origin$Origin <- factor(origin$Origin, levels=c("nLung",'tLung','nLN','mLN'))
for(c in levels(origin$Cluster)) origin$Number[origin$Cluster==c] <- origin$Number[origin$Cluster==c]/sum(origin$Number[origin$Cluster==c])*100

p2 <- ggplot(origin, aes(x=Cluster, y=Number, fill=Origin)) +
    geom_bar(stat="identity", width=0.8) + 
    scale_fill_manual(values=c("gold1","firebrick","darkseagreen2","green4")) +
    scale_y_continuous(position = "right") +
    labs(x=NULL, y="Percentage", size=.5)+
    theme(panel.background=element_blank())

pdf('barplot_kim_llb_b25_stack.pdf', width=6.5, height=2.5)
p2
dev.off()

### Documentation ####

library(Seurat)
?Seurat ## v3.1.0
rm(list=ls())

# load('kim_umi_llb_seurat.RData')
load('kim_umi_ll_seurat.RData')

### Check Cluster 14 of CD4 for IL17A and IL13 #### 

rm(list=ls())
# load('kim_umi_ll_seurat.RData')
load('kim_umi_llb_seurat.RData')

il17a <- kp[['RNA']]['IL17A', kp$seurat_clusters==14]
il13 <- kp[['RNA']]['IL13', kp$seurat_clusters==14]

# plot(density(il17a %>% as.numeric))
# plot(density(il13 %>% as.numeric))
summary(il17a %>% as.numeric)
summary(il13 %>% as.numeric)
sum(il17a>0)
sum(il13>0)
sum(kp$seurat_clusters==14)

ncell <- table(kp$seurat_clusters) %>% as.data.frame
names(ncell) <- c("cluster", "#cell")
write.csv(ncell, file="cluster_number_b.csv", row.names=F)


### write cluster genes ####
ll <- read.csv('top30_ll.csv')
llb <- read.csv('top30_llb.csv')
library(xlsx)
ll <- ll[c(8,3:7)]
names(ll) <- c("Gene","Log FoldChange","Propotion positive in cluster",
               "Proportion positive outside cluster","FDR-corrected P-value","Cluster")
write.xlsx(ll, file="top30_byCluster.xlsx", sheetName = "T cells", row.names = F)
llb <- llb[c(8,3:7)]
names(llb) <- c("Gene","Log FoldChange","Propotion positive in cluster",
               "Proportion positive outside cluster","FDR-corrected P-value","Cluster")
write.xlsx(llb, file="top30_byCluster.xlsx", sheetName = "B cells", row.names=F, append=T)





##### TCR overlap anlaysis (For Connelli, 21/5/13) ####
setwd('~/Documents/Projects/cc/kelli/')
library(dplyr)
library(stringr)
rm(list=ls())
clones_raw <- xlsx::read.xlsx('frequency of clones of total.xlsx', sheetIndex=1)
clones_raw$NA..3 %>% head
clones_sum <- clones_raw[1:4, 5:7]
clones <- clones_raw[1:4]; names(clones) <- c('Group','Clonetype','Freq','TRAB')
clones$Clonetype <- clones$Clonetype %>% as.character %>% gsub('clonotype','',.) %>% as.numeric
clones$TRA <- clones_raw$NA..3 %>% as.character %>% str_extract(pattern = "TRA:[A-Z]+_") %>% gsub('TRA:|_','',.)
clones$TRB <- clones_raw$NA..3 %>% as.character %>% str_extract(pattern = "TRB:[A-Z]+") %>% gsub('TRB:','',.)
clones$TRA[is.na(clones$TRA)] <- ""
clones$TRB[is.na(clones$TRB)] <- ""
save(clones, clones_sum, file="TCRclones.RData")

### simpson's diversity index
types <- clones$TRAB %>% unique
byclone <- data.frame(type=types)
byclone$wk8_tumor <- clones$Freq[clones$Group=="8wk_tumor"][match(types, clones$TRAB[clones$Group=="8wk_tumor"])]
byclone$wk8_tumor[is.na(byclone$wk8_tumor)] <- 0
byclone$wk8_LN <- clones$Freq[clones$Group=="8wk_LN"][match(types, clones$TRAB[clones$Group=="8wk_LN"])]
byclone$wk8_LN[is.na(byclone$wk8_LN)] <- 0
byclone$wk17_tumor <- clones$Freq[clones$Group=="17wk_tumor"][match(types, clones$TRAB[clones$Group=="17wk_tumor"])]
byclone$wk17_tumor[is.na(byclone$wk17_tumor)] <- 0
byclone$wk17_mLN <- clones$Freq[clones$Group=="17wk_mLN"][match(types, clones$TRAB[clones$Group=="17wk_mLN"])]
byclone$wk17_mLN[is.na(byclone$wk17_mLN)] <- 0
byclone$p8tm <- byclone$wk8_tumor/sum(byclone$wk8_tumor)
byclone$p8ln <- byclone$wk8_LN/sum(byclone$wk8_LN)
byclone$p17tm <- byclone$wk17_tumor/sum(byclone$wk17_tumor)
byclone$p17ln <- byclone$wk17_mLN/sum(byclone$wk17_mLN)
si_8wk_tm <- byclone$p8tm %>% .^2 %>% sum
si_8wk_ln <- byclone$p8ln %>% .^2 %>% sum
si_17wk_tm <- byclone$p17tm %>% .^2 %>% sum
si_17wk_ln <- byclone$p17ln %>% .^2 %>% sum
mh_tm = 2*sum(byclone$p8tm*byclone$p17tm)/(si_8wk_tm + si_17wk_tm)
mh_ln = 2*sum(byclone$p8ln*byclone$p17ln)/(si_8wk_ln + si_17wk_ln)
mh_8wk = 2*sum(byclone$p8tm*byclone$p8ln)/(si_8wk_tm + si_8wk_ln)
mh_17wk = 2*sum(byclone$p17tm*byclone$p17ln)/(si_17wk_tm + si_17wk_ln)
mh_8tm_17ln = 2*sum(byclone$p8tm*byclone$p17ln)/(si_8wk_tm + si_17wk_ln)
mh_8ln_17tm = 2*sum(byclone$p17tm*byclone$p8ln)/(si_17wk_tm + si_8wk_ln)


write.csv(clones_sum, file="clones_stats.csv", row.names=F)



###### THE END #######
