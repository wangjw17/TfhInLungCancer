#### Tfh/B cell proportion analysis 200510
setwd('~/Documents/Projects/cc/tfh_fraction/')
library(dplyr)

### check TCGA downloaded data ====
meta <- read.table('tcga_luad_top100/MANIFEST.txt', header=T)
meta$filename <- as.character(meta$filename)
meta$filename <- gsub(".gz","",meta$filename)
meta$fname <- meta$filename %>% as.character %>% strsplit(split="/") %>% sapply("[",2)
meta$norm <- NA
meta$norm[grep('htseq',meta$fname)] <- "Counts"
meta$norm[grep('FPKM', meta$fname)] <- "FPKM"
meta$norm[grep('FPKM-UQ', meta$fname)] <- "FPKM-UQ"
samples <- dir(".")
samples <- samples[-grep("MANIFEST",samples)]
sams <- meta$fname %>% strsplit(split="-") %>% sapply("[",5) #%>% gsub(".txt","",.)
for (i in 1:100){
  mat <- read.table(meta$filename[i], header=T)
  names(mat) <- c("Geneid",sams[i])
  if (exists("matall")) matall <- left_join(matall, mat, by="Geneid") else matall <- mat
}
dim(matall)
View(matall)
meta$filename <- as.character(meta$filename)
abline(0,1)
## BG: FPKM = [RMg * 109 ] / [RM75 * L]; scale difference with FPKM
plot(rowMeans(matall[,1+which(meta$norm=="FPKM")]), rowMeans(matall[,1+which(meta$norm=="FPKM-UQ")]), 
     log="xy", xlab="FPKM", ylab="FPKM-UQ")
plot(matall$`710af9d30710.FPKM.txt`,matall$`5adb24f71754.FPKM.txt`, log="xy")
matall1 <- matall[,c(1,1+which(meta$norm!="Counts"))]
save(matall1, meta, file="FPKM_top100.RData")
hist(rowMeans(matall1[,-1]) %>% log)
matall2 <- matall1[rowMeans(matall1[,-1])>1,]
dim(matall2)
matall2$Geneid <- matall2$Geneid %>% strsplit(split="\\.") %>% sapply("[",1)
write.table(matall2, file="mixture_luad.txt", sep="\t", quote=F, row.names = F)
matall2$Geneid <- probes$external_gene_name[match(matall2$Geneid,probes$ensembl_gene_id)]
write.table(matall2, file="mixture_luad_h.txt", sep="\t", quote=F, row.names = F)


### download signature matrix ====
## LM22 (CIBERSORT has it)
## Immunostates
imms <- xlsx::read.xlsx('~/Downloads/41467_2018_7242_MOESM5_ESM.xlsx', sheetIndex = 1)
header <- imms[1,] %>% as.matrix %>% as.character 
imms <- imms[-1,]; names(imms) <- header
View(imms)
load('../../PTSD/data/datProbes.RData')
ens <- probes$ensembl_gene_id[match(imms$Gene,probes$external_gene_name)]
imms$Gene <- ens
imms <- imms[!is.na(imms$Gene),]
write.table(imms, file="sigmat_imms.txt", sep="\t", quote=F, row.names=F)

### read/plot result files ====
library(ggplot2)
library(tidyr)
library(reshape2)
res <- read.csv('~/Downloads/CIBERSORT.Output_Job6.csv')
boxplot(res[,2:21])
meds <- apply(res[,2:21],2,median)
cells <- names(res)[2:21]
resp <- melt(res[,1:21])
resp$variable <- factor(resp$variable, levels=cells[order(meds, decreasing = T)])
pdf('fractions_immunos.pdf', width=12)
p <- ggplot(resp, aes(x=variable, y=value)) + 
  geom_boxplot() + 
  labs(title = "Cell type fractions based on ImmunoStates") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 
print(p)
dev.off()

# res <- read.csv('~/Downloads/CIBERSORT.Output_Job7.csv')
res <- read.csv('~/Downloads/CIBERSORT.Output_Job8.csv')
boxplot(res[,2:23])
meds <- apply(res[,2:23],2,median)
cells <- names(res)[2:23]
resp <- melt(res[,1:23])
resp$variable <- factor(resp$variable, levels=cells[order(meds, decreasing = T)])
pdf('fractions_lm22.pdf', width=12)
p <- ggplot(resp, aes(x=variable, y=value)) + 
  geom_boxplot() + 
  labs(title = "Cell type fractions based on LM22, Primary Tumor") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 
print(p)
dev.off()

### plot of combined cell types ====
# 1. B cells: B cells naïve; B cells memory; Plasma cells
# 2. CD8 T cells
# 3. CD4 naïve T cells
# 4. CD4 memory T cells: T cells CD4 memory resting; T cells CD4 memory activated
# 5. T follicular helper cells
# 6. Regulatory T cells
# 7. γδ T cells: T cell gamma delta
# 8. NK cells: NK cells resting; NK cells activated
# 9. Monocytes
# 10. Macrophages: Macrophages M0; M1; M2
# 11. Dendritic cells: DC resting; DC activated
# 12. Mast cells: Mast cells resting; Mast cells activated
# 13. Eosinophils
# 14. Neutrophils
library(reshape2)
library(ggplot2)
# res <- read.csv('~/Downloads/CIBERSORT.Output_Job8.csv')
res <- read.csv('CIBERSORT.Output_Job8.csv')
res2 <- res[c(5:6,9:11,14,22:23)]
names(res2)[4] <- "T.cells.regulatory"
res2$B.cells <- rowSums(res[c(2:4)])
res2$T.cells.CD4.memory <- rowSums(res[c(7:8)])
res2$NK.cells <- rowSums(res[c(12:13)])
res2$Macrophage <- rowSums(res[c(15:17)])
res2$Dendritic.cells <- rowSums(res[18:19])
res2$Mast.cells <- rowSums(res[20:21])
resp <- melt(cbind(res[1],res2))
meds <- apply(res2,2,median)
resp$variable <- factor(resp$variable, levels=names(res2)[order(meds, decreasing = T)])
pdf('fractions_lm22_ptumor_comb.pdf', width=12)
p <- ggplot(resp, aes(x=variable, y=value)) + 
  geom_boxplot() +
  labs(title = "Cell type fractions based on LM22, Primary Tumor") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 
print(p)
dev.off()

pdf('fractions_lm22_ptumor_bar.pdf', width=12)
resb <- colMeans(res2) %>% as.data.frame
names(resb)[1] <- "fraction"; resb$cell <- rownames(resb)
resb$cell <- factor(resb$cell, levels=resb$cell[order(resb$fraction, decreasing = T)])
p <- ggplot(resb, aes(x=cell,y=fraction)) + 
  geom_bar(stat="identity") +
  labs(title = "Cell type fractions based on LM22, Primary Tumor") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 
print(p)
dev.off()





### process 1782 files ====
tb <- read.table('gdc_manifest_20200520_043845.txt', header=T)
tb$filename %>% grep("FPKM.txt.gz",.) %>% length
tb$filename %>% grep("FPKM-UQ.txt.gz",.) %>% length
tb$filename %>% grep("htseq.counts.gz",.) %>% length

## there're (594 Files/585 Samples/515 Cases), each with FPKM, FPKM-UQ, htseq.counts; use FPKM only
tb1 <- read.table('gdc_sample_sheet.2020-05-19.tsv', header=T, sep="\t")
table(tb1$Case.ID) %>% length; table(tb1$Case.ID) %>% max
table(tb1$Sample.ID) %>% length; table(tb1$Sample.ID) %>% max
id.sel <- intersect(which(tb1$Sample.Type == "Primary Tumor"), grep('FPKM.txt',tb1$File.Name))
View(tb1[tb1$Case.ID=="TCGA-44-2668",])
## use only one sample from each case
tb2 <- tb1[id.sel,]
id.sel2 <- sapply(tb2$Case.ID %>% unique, function(a) grep(a, tb2$Case.ID)[1])
tb3 <- tb2[id.sel2,]
##> filters: 0) FPKM --> 594 --> 1) Primary Tumor --> 533 --> 2) First sample per CaseID --> 513
fnames <- paste0(tb3$File.ID, "/", tb3$File.Name)
for (fname in fnames) system(paste0('mv gdc_download_20200520_040142.574884/',fname, ' gdc_luad_ptumor_fpkm/'))
system("gunzip gdc_luad_ptumor_fpkm/*")

### concatenate 513 files ====
fnames <- dir("gdc_luad_ptumor_fpkm/")
fshort <- fnames %>% strsplit(split="-") %>% sapply("[",1) ## it's unique
for (i in 1:length(fnames)){
  f <- fnames[i]; ff <- fshort[i]
  mat <- read.table(paste0('gdc_luad_ptumor_fpkm/',f), header=T)
  names(mat) <- c("Geneid",ff)
  if (exists("matall")) matall <- left_join(matall, mat, by="Geneid") else matall <- mat
}
dim(matall)
View(matall)
hist(rowMeans(matall[,-1]) %>% log1p)
matall1 <- matall[rowMeans(matall[,-1])>1,] ## filter the lowest genes
dim(matall1)
matall1$Geneid <- matall1$Geneid %>% as.character %>% strsplit(split="\\.") %>% sapply("[",1)
load('../../../PTSD/data/datProbes.RData')
matall1$Geneid <- probes$external_gene_name[match(matall1$Geneid,probes$ensembl_gene_id)]
matall2 <- matall1[!is.na(matall1$Geneid),] ## filter the unmatched genes
write.table(matall2, file="gdc_luad_ptumor_fpkm/mixture_luad_tpkm.txt", sep="\t", quote=F, row.names = F)

## 2020/10/27 save original and processed data
matall$Geneid <- matall$Geneid %>% as.character %>% strsplit(split="\\.") %>% sapply("[",1)
matall$Genename <- probes$external_gene_name[match(matall$Geneid,probes$ensembl_gene_id)]
save(matall, matall2, file="mixture_luad_fpkm.RData")

### downlaod MC38 data ====
setwd("../mc38/")
tb1 <- read.table('GSE124691_Genes.tsv', header=F, sep="\t")


### published pipeline ====
BiocManager::install("limma")
BiocManager::install("grimbough/biomaRt")
install.packages("robustSingleCell")

### use Seurat to read in data (TIL exp1; dLN exp1) ====
setwd("~/Documents/Projects/cc/mc38/")
library(Seurat)

## read in data
# rawdata <- Read10X(data.dir = "GSE124691_RAW/GSM3543446_Experiment_1_TILs_CD44pos_PDpos/")
rawdata <- Read10X(data.dir = "GSE124691_RAW/GSM3543443_Experiment_1_dLN_GP66pos_merged/") ##20200719
kp <- CreateSeuratObject(counts = rawdata, project = "GSM3543446_Exp1", min.cells = 3, min.features = 200)
summary(kp)
kp[["percent.mt"]] <- PercentageFeatureSet(kp, pattern = "^mt-")
VlnPlot(kp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# plot1 <- FeatureScatter(kp, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(kp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# CombinePlots(plots = list(plot1, plot2))
kp <- subset(kp, subset = nCount_RNA < 20000 & nFeature_RNA < 3000 & percent.mt < 5)

## normalization
kp <- NormalizeData(kp, normalization.method = "LogNormalize", scale.factor = 10000)
kp <- FindVariableFeatures(kp, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(kp), 10)
plot1 <- VariableFeaturePlot(kp)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
CombinePlots(plots = list(plot1, plot2))

## dimension reduction
all.genes <- rownames(kp)
kp <- ScaleData(kp, features = all.genes)
kp <- RunPCA(kp, features = VariableFeatures(kp))
print(kp[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(kp, dims = 1:5, reduction = "pca")
DimHeatmap(kp, dims = 1:12, cells = 500, balanced = TRUE)
DimPlot(kp, reduction="pca")
# ## determine dimensionality
# kp <- JackStraw(kp, num.replicate = 100) 
# kp <- ScoreJackStraw(kp, dims = 1:20)
# JackStrawPlot(kp, dims = 1:20)
# ElbowPlot(kp)
## clustering cells
kp <- FindNeighbors(kp, dims = 1:20)
kp <- FindClusters(kp, resolution = 0.5)
head(Idents(kp), 5)
## UMPA
kp <- RunUMAP(kp, dims = 1:20)
DimPlot(kp, reduction = "umap")

## plot by markers
kp.markers <- FindAllMarkers(kp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top20 <- kp.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
write.csv(top20, file="results/top_genes_dln.csv", row.names=F)
# top50 <- kp.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
# saveRDS(kp.markers, file = "Luc14_markers_f3.rds")
# DoHeatmap(kp, features = top20$gene) + NoLegend()
pdf('results/marker_expression_dln.pdf', width=10, height=10)
DimPlot(kp, reduction = "umap")
VlnPlot(kp, features = c("Ccr7", "Cxcr3", "Cxcr5","Lef1","Tcf7","Bcl6"))
FeaturePlot(kp, features = c("Ccr7", "Cxcr3", "Cxcr5","Lef1","Tcf7","Bcl6"))
dev.off()


### Kumar and Perou data ====
## Kumar Data GSE121861
library(ggplot2)
setwd('../kumar')
rm(list=ls())
rows <- read.csv('GSE121861_syngeneic_row_data.csv')
cols <- read.csv('GSE121861_syngeneic_column_data.csv')
expmat <- read.csv('GSE121861_syngeneic_expression.csv')
sel <- cols$Protocol=="fresh CD45pos"
expmat <- expmat[sel,] %>% t()
cols <- cols[sel,]
gene.sel <- rowMeans(expmat>0, na.rm = T)>.01
rows <- rows[gene.sel,]
expmat <- expmat[gene.sel,]
rm.na <- colSums(is.na(expmat))>0
expmat <- expmat[,!rm.na]
cols <- cols[!rm.na,]

## dimension reduction
pca.tpm <- prcomp(t(expmat), scale. = F, center = T)
pcatpm <- pca.tpm$x
plot(pcatpm[,1:2])
tsntpm <- tsne::tsne(pcatpm[,1:20], k = 2)
plot(tsntpm, pch=16)
km <- kmeans(tsntpm, 10)
tsntpm <- as.data.frame(tsntpm)
names(tsntpm) <- c("tSNE1","tSNE2")
ggplot(tsntpm, aes(tSNE1,tSNE2,col=factor(km$cluster)))+
  geom_point()
umptpm <- umap::umap(pcatpm[,1:20])$layout
umptpm <- as.data.frame(umptpm)
names(umptpm) <- c("UMAP1","UMAP2")
id.um <- kmeans(umptpm, 10)$cluster
ggplot(umptpm, aes(UMAP1,UMAP2,col=factor(id.um)))+
  geom_point()

library(Seurat)
rownames(expmat) <- rows$GeneSymbol
kp <- CreateSeuratObject(counts = expmat, project = "KumarDataset", min.cells = 3, min.features = 200)
summary(kp)
kp[["percent.mt"]] <- PercentageFeatureSet(kp, pattern = "^mt-")
VlnPlot(kp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# plot1 <- FeatureScatter(kp, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(kp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# CombinePlots(plots = list(plot1, plot2))
kp <- subset(kp, subset = nCount_RNA < 15000 & nFeature_RNA < 5000 & percent.mt < 10)


## Perou data
setwd('../perou/')
rm(list=ls())
expmat <- read.table('GSE124821_In_Vivo_Normalized_Matrix.txt', sep="\t")
meta <- expmat[1:23,] %>% t()
colnames(meta) <- meta[1,]
meta <- meta[-1,]
datmat <- expmat[27:nrow(expmat),-1]
cols<- expmat[26,-1]
rows <- expmat[27:nrow(expmat),1]
dat <- apply(datmat,2,function(x) as.numeric(as.matrix(x)))
colnames(dat) <- cols
rownames(dat) <- rows
dat <- dat[,-265]

library(Seurat)
kp <- CreateSeuratObject(counts = dat, project = "PerouDataset", min.cells = 3, min.features = 200)
summary(kp)
kp[["percent.mt"]] <- PercentageFeatureSet(kp, pattern = "^mt-")
VlnPlot(kp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kp <- subset(kp, subset = nFeature_RNA > 20000)


## Common for Kumar and Perou
## normalization
kp <- NormalizeData(kp, normalization.method = "LogNormalize", scale.factor = 10000)
kp <- FindVariableFeatures(kp, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(kp), 10)
plot1 <- VariableFeaturePlot(kp)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
CombinePlots(plots = list(plot1, plot2))

## dimension reduction
all.genes <- rownames(kp)
kp <- ScaleData(kp, features = all.genes)
kp <- RunPCA(kp, features = VariableFeatures(kp))
print(kp[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(kp, dims = 1:5, reduction = "pca")
DimHeatmap(kp, dims = 1:12, cells = 500, balanced = TRUE)
DimPlot(kp, reduction="pca")
# ## determine dimensionality
# kp <- JackStraw(kp, num.replicate = 100) 
# kp <- ScoreJackStraw(kp, dims = 1:20)
# JackStrawPlot(kp, dims = 1:20)
# ElbowPlot(kp)
## clustering cells
kp <- FindNeighbors(kp, dims = 1:20)
kp <- FindClusters(kp, resolution = 0.5)
head(Idents(kp), 5)
## UMPA
kp <- RunUMAP(kp, dims = 1:20)
DimPlot(kp, reduction = "umap")

## plot by markers
kp.markers <- FindAllMarkers(kp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top20 <- kp.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
write.csv(top20, file="top_genes_perou.csv", row.names=F)
# top50 <- kp.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
# saveRDS(kp.markers, file = "Luc14_markers_f3.rds")
# DoHeatmap(kp, features = top20$gene) + NoLegend()
pdf('marker_expression_perou.pdf', width=10, height=10)
DimPlot(kp, reduction = "umap")
VlnPlot(kp, features = c("Ccr7", "Cxcr3", "Cxcr5","Lef1","Tcf7","Bcl6"))
FeaturePlot(kp, features = c("Ccr7", "Cxcr3", "Cxcr5","Lef1","Tcf7","Bcl6"))
dev.off()


### replot by Can ====
rm(list=ls())

library(ggplot2)
library(reshape2)
library(dplyr)
library(hrbrthemes)

res <- read.csv('CIBERSORT.Output_Job8.csv')
res2 <- res[c(5:6,9:11,14,22:23)]
names(res2)[1] <- "CD8.T.cells"
names(res2)[2] <- "CD4.naive.T.cells"
names(res2)[3] <- "T.follicular.helper.cells"
names(res2)[4] <- "Regulatory.T.cells"
names(res2)[5] <- "γδ.T.cells"
res2$B.cell.lineage <- rowSums(res[c(2:4)])
res2$CD4.memory.T.cells <- rowSums(res[c(7:8)])
res2$NK.cells <- rowSums(res[c(12:13)])
res2$Macrophage <- rowSums(res[c(15:17)])
res2$Dendritic.cells <- rowSums(res[18:19])
res2$Mast.cells <- rowSums(res[20:21])
resp <- melt(cbind(res[1],res2))
meds <- apply(res2,2,median)
resp$variable <- factor(resp$variable, levels=names(res2)[order(meds, decreasing = T)])

p <- ggplot(resp1, aes(x=variable, y=value, fill=variable)) + 
  geom_violin(scale=3) +
  # geom_boxplot(width=0.1) +
  theme_classic() +
  labs(title = "Cell type fractions based on LM22, Primary Tumor") +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12, color = "black"),
        axis.text.y = element_text(angle = 0, size=12, color = "black")) +
  scale_fill_manual(values=c("green4", "gray96", "gray96", "gray96", "gray96","gray96", "gray96",
                             "firebrick", "gray96", "gray96","gray96", "gray96","gray96", "gray96"))
p

### plot GCB/CD4T subspecies ====
pdf('fractions_lm22_ptumor_comb_violin.pdf', width=12, height=8)
p
dev.off()

resp1 <- res[c('Input.Sample','B.cells.naive','B.cells.memory','Plasma.cells')] %>% melt
p <- ggplot(resp1, aes(x=variable, y=value, fill=variable)) + 
  geom_violin(scale=3) +
  # geom_boxplot(width=0.1) +
  theme_classic() +
  labs(title = "Cell type fractions based on LM22, Primary Tumor") +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12, color = "black"),
        axis.text.y = element_text(angle = 0, size=12, color = "black"))
p
pdf('fractions_lm22_ptumor_GCB_violin.pdf', width=6, height=6)
p
dev.off()

resp1 <- res[c('Input.Sample','T.cells.CD4.naive','T.cells.follicular.helper','T.cells.regulatory..Tregs.',
               'T.cells.CD4.memory.resting','T.cells.CD4.memory.activated')] %>% melt
resp1$variable <- gsub(".resting|.activated", "", resp1$variable)
p <- ggplot(resp1, aes(x=variable, y=value, fill=variable)) + 
  geom_violin(scale=3) +
  # geom_boxplot(width=0.1) +
  theme_classic() +
  labs(title = "Cell type fractions based on LM22, Primary Tumor") +
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12, color = "black"),
        axis.text.y = element_text(angle = 0, size=12, color = "black"))
p
pdf('fractions_lm22_ptumor_CD4T_violin.pdf', width=6, height=6)
p
dev.off()



### 1. TCGA LUAD: A=longer survival; B=shorter survival =====
#   1) Compare TFH/GC B cell fraction in A vs. B
#   2) Compare IL21 level in A vs. B
tb <- read.table('gdc_manifest_20200520_043845.txt', header=T)
tb$filename %>% grep("FPKM.txt.gz",.) %>% length
tb$filename %>% grep("FPKM-UQ.txt.gz",.) %>% length
tb$filename %>% grep("htseq.counts.gz",.) %>% length

## there're (594 Files/585 Samples/515 Cases), each with FPKM, FPKM-UQ, htseq.counts; use FPKM only
tb1 <- read.table('gdc_sample_sheet.2020-05-19.tsv', header=T, sep="\t")
table(tb1$Case.ID) %>% length; table(tb1$Case.ID) %>% max
table(tb1$Sample.ID) %>% length; table(tb1$Sample.ID) %>% max
id.sel <- intersect(which(tb1$Sample.Type == "Primary Tumor"), grep('FPKM.txt',tb1$File.Name))
View(tb1[tb1$Case.ID=="TCGA-44-2668",])
## use only one sample from each case
tb2 <- tb1[id.sel,]
id.sel2 <- sapply(tb2$Case.ID %>% unique, function(a) grep(a, tb2$Case.ID)[1])
tb3 <- tb2[id.sel2,]

##! only 267 total annotations.txt (FPKM) in downloaded folder; only 89 are FPKM
##! only 70 have matched ids and survival info

## using aliquot and clinical tsv
manif <- read.table('gdc_download_20200520_040142.574884/MANIFEST.txt', header=T)
cbx <- read.csv('CIBERSORT.Output_Job8.csv')
aliq <- read.table('biospecimen.cart.2020-05-19/aliquot.tsv', sep="\t", header=T, quote="\"")[,1:11]
clic <- read.table('clinical.cart.2020-05-19/clinical.tsv', sep="\t", quote="\"", header=T)
clic <- apply(clic, 2, function(x) gsub("'--",NA,x)) %>% as.data.frame
clic <- clic[,colSums(!is.na(clic))>0.25*nrow(clic)]
clic$survival <- clic$days_to_death %>% as.character %>% as.numeric
clic$survival <- clic$survival/365
sum(!is.na(clic$survival)) ## 366 samples with survival information and case_id/case_submitter_id
clic$case_id %in% aliq$case_id
samples <- dir("gdc_download_20200520_040142.574884/")

## CBx ID = tb3$File.name <=> tb3$Case.ID = clic$Case_submitter_id <=> clic$case_id
tsam <- cbx$Input.Sample[1] %>% as.character
csam <- clic$case_id[1] %>% as.character %>% strsplit(split="-") %>% sapply("[",1)
files <- dir("gdc_download_20200520_040142.574884/")
fsam <- files[1] %>% as.character
##> manifest files (tb, manif) only have folder + filename info
tb3$Sample_CBx <- tb3$File.Name %>% as.character %>% strsplit(split="-") %>% sapply("[",1)
tb3$survival <- clic$survival[match(tb3$Case.ID, clic$case_submitter_id)]
## 183 samples out of 513 have survival information
write.csv(tb3, file = "gdc_sample_summary_1027.csv", row.names=F)

### plot by higher or lower survival time ====
library(ggplot2)
library(reshape2)
library(dplyr)
library(hrbrthemes)

res <- read.csv('CIBERSORT.Output_Job8.csv')
res2 <- res[c(5:6,9:11,14,22:23)]
names(res2)[1] <- "CD8.T.cells"
names(res2)[2] <- "CD4.naive.T.cells"
names(res2)[3] <- "T.follicular.helper.cells"
names(res2)[4] <- "Regulatory.T.cells"
names(res2)[5] <- "γδ.T.cells"
res2$B.cell.lineage <- rowSums(res[c(2:4)])
res2$CD4.memory.T.cells <- rowSums(res[c(7:8)])
res2$NK.cells <- rowSums(res[c(12:13)])
res2$Macrophage <- rowSums(res[c(15:17)])
res2$Dendritic.cells <- rowSums(res[18:19])
res2$Mast.cells <- rowSums(res[20:21])

cbx <- read.csv('CIBERSORT.Output_Job8.csv', header=T)
cbx$survival <- tb3$survival[match(cbx$Input.Sample, tb3$Sample_CBx)]
sel.surv <- !is.na(cbx$survival)
res2 <- res2[sel.surv,]
surv_med <- median(cbx$survival, na.rm=T)
surv_top <- quantile(cbx$survival, 0.75, na.rm=T)
surv_bot <- quantile(cbx$survival, 0.25, na.rm=T)
cbx$surv_time <- NA
cbx$surv_time[cbx$survival>=surv_top] <- "long"
cbx$surv_time[cbx$survival<surv_bot] <- "short"
# cbx <- cbx[sel.surv,]

resp <- melt(cbind(res[sel.surv,1],res2))
names(resp)[1] <- "Input.Sample"
meds <- apply(res2,2,median)
resp$variable <- factor(resp$variable, levels=names(res2)[order(meds, decreasing = T)])
resp$survival <- cbx$surv_time[match(resp$Input.Sample,cbx$Input.Sample)]
resp <- resp[!is.na(resp$survival),]


p <- ggplot(resp, aes(x=variable, y=value, fill=survival)) + 
  geom_boxplot() + 
  theme_classic() +
  labs(title = "Cell type fractions based on LM22, Primary Tumor") +
  xlab("") +
  ylab("")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12, color = "black"),
        axis.text.y = element_text(angle = 0, size=12, color = "black")) +
  scale_fill_manual(values=c("green4", "gray96", "gray96", "gray96", "gray96","gray96", "gray96", 
                             "firebrick", "gray96", "gray96","gray96", "gray96","gray96", "gray96"))
p
pdf('fractions_lm22_ptumor_comb_surv_quartile.pdf', height=8, width=12)
print(p)
dev.off()
write.csv(cbx, file="CIBERSORT.Output_Job8_surv.csv")


### plot IL21 and CXCL13 by survival time ====
rm(list=ls())
tb3 <- read.csv('gdc_sample_summary_1027.csv')
load('mixture_luad_fpkm.RData')
# matall <- read.table('mixture_luad_tpkm.txt', header=T)
# load('../../../PTSD/data/datProbes.RData')
# ensg.il21 <- datProbes$ensembl_gene_id[datProbes$external_gene_id=="IL21"]
# ensg.cxcl <- datProbes$ensembl_gene_id[datProbes$external_gene_id=="CXCL13"]
## IL21 not in the matrix
df <- data.frame(Sample=colnames(matall)[-c(1,515)] %>% gsub("X","",.), 
                 IL21=matall[which(matall$Genename=="IL21"),-c(1,515)] %>% as.numeric,
                 CXCL13=matall[which(matall$Genename=="CXCL13"),-c(1,515)] %>% as.numeric)
df$survival <- tb3$survival[match(df$Sample, tb3$Sample_CBx)]
df <- df[!is.na(df$survival),]
df$surv_half <- as.numeric(df$survival >= median(df$survival)) %>% factor
df$surv_quar <- as.numeric(df$survival >= quantile(df$survival, 0.75))
df$surv_quar[df$survival <= quantile(df$survival, 0.25)] <- -1
df$surv_quar <- df$surv_quar %>% factor
df$surv_half <- as.numeric(df$surv_half) %>% as.factor
df$surv_four <- 2
df$surv_four[df$surv_quar==-1] <- 1
df$surv_four[df$surv_quar==1] <- 4
df$surv_four[df$surv_half==2 & df$surv_four==2] <- 3
df$surv_four <- as.factor(df$surv_four)
boxplot(df$CXCL13 %>% log1p ~df$surv_half)
boxplot(df$CXCL13 %>% log1p ~df$surv_quar)


library(ggplot2)
pdf('expression_ptumor_IL21_CXCL13_violin.pdf')
ggplot(df, aes(x=surv_half, y=log1p(IL21), fill=surv_half)) +
  geom_violin() +
  xlab('quarter, survival time') + ylab('logFPKM, IL21') +
  labs(title="Expression of IL21 in survival time-stratified groups") +
  scale_fill_manual(values = c("grey40","grey80")) +
  theme_classic()
ggplot(df, aes(x=surv_quar, y=log1p(IL21), fill=surv_quar)) +
  geom_violin()+
  xlab('quarter, survival time') + ylab('logFPKM, IL21') +
  labs(title="Expression of IL21 in survival time-stratified groups") +
  scale_fill_manual(values = c("grey40","grey60", "grey80")) +
  theme_classic()
ggplot(df, aes(x=surv_four, y=log1p(IL21), fill=surv_four)) +
  geom_violin()+
  xlab('quarter, survival time') + ylab('logFPKM, IL21') +
  labs(title="Expression of IL21 in survival time-stratified groups") +
  scale_fill_manual(values = c("grey40","grey50", "grey70", "grey80")) +
  theme_classic()
ggplot(df[df$surv_four%in%c(1,4),], aes(x=surv_four, y=log1p(IL21), fill=surv_four)) +
  geom_violin()+
  xlab('quarter, survival time') + ylab('logFPKM, IL21') +
  labs(title="Expression of IL21 in survival time-stratified groups") +
  scale_fill_manual(values = c("grey40","grey80")) +
  theme_classic()
ggplot(df, aes(x=surv_half, y=log1p(CXCL13),fill=surv_half)) +
  geom_violin()+
  xlab('quarter, survival time') + ylab('logFPKM, CXCL13') +
  labs(title="Expression of CXCL13 in survival time-stratified groups") +
  scale_fill_manual(values = c("grey40","grey80")) +
  theme_classic()
ggplot(df, aes(x=surv_quar, y=log1p(CXCL13), fill=surv_quar)) +
  geom_violin()+
  xlab('quarter, survival time') + ylab('logFPKM, CXCL13') +
  labs(title="Expression of CXCL13 in survival time-stratified groups") +
  scale_fill_manual(values = c("grey40","grey60", "grey80")) +
  theme_classic()
ggplot(df, aes(x=surv_four, y=log1p(CXCL13), fill=surv_four)) +
  geom_violin()+
  xlab('quarter, survival time') + ylab('logFPKM, CXCL13') +
  labs(title="Expression of CXCL13 in survival time-stratified groups") +
  scale_fill_manual(values = c("grey40","grey50", "grey70", "grey80")) +
  theme_classic()
ggplot(df[df$surv_four%in%c(1,4),], aes(x=surv_four, y=log1p(CXCL13), fill=surv_four)) +
  geom_violin()+
  xlab('quarter, survival time') + ylab('logFPKM, CXCL13') +
  labs(title="Expression of CXCL13 in survival time-stratified groups") +
  scale_fill_manual(values = c("grey40","grey80")) +
  theme_classic()
dev.off()

library(ggpubr)
ggplot(df[df$surv_four%in%c(1,4),], aes(x=surv_four, y=log1p(CXCL13), fill=surv_four)) +
  geom_point(position="jitter")+
  # geom_boxplot() +
  xlab('quarter, survival time') + ylab('logFPKM, CXCL13') +
  labs(title="Expression of CXCL13 in survival time-stratified groups") +
  scale_fill_manual(values = c("grey40","grey80")) +
  stat_compare_means() +
  theme_classic()
ggplot(df[df$surv_four%in%c(1,4),], aes(x=surv_four, y=log1p(IL21), fill=surv_four)) +
  geom_violin()+
  xlab('quarter, survival time') + ylab('logFPKM, IL21') +
  labs(title="Expression of IL21 in survival time-stratified groups") +
  scale_fill_manual(values = c("grey40","grey80")) +
  stat_compare_means() +
  theme_classic()

wilcox.test(log1p(df$CXCL13[df$survival<quantile(df$survival,1/3)]), 
            log1p(df$CXCL13[df$survival>quantile(df$survival,2/3)]))
wilcox.test(log1p(df$IL21[df$survival<quantile(df$survival,1/3)]), 
            log1p(df$IL21[df$survival>quantile(df$survival,2/3)]))

write.csv(df, file="expression_ptumor_survival_IL21_CXCL13.csv", row.names=F)

# library(reshape2)
# df2 <- melt(df[df$surv_four%in%c(1,4), c("surv_four","IL21","CXCL13")])
# ggplot(df2, aes(x=variable, y=log1p(value), fill=surv_four)) +
#   geom_violin()+
#   xlab('marker x survival time') + ylab('logFPKM') +
#   labs(title="Expression of IL21 and CXCL13 in survival time-stratified groups") +
#   scale_fill_manual(values = c("grey40","grey80")) +
#   theme_classic()


### 2. Is there correlation of IL21 with survival? Using TCGA-LUAD cohort (X) ====
### 3. Analyze scRNA seq data from Kim et.al. nature communication 2020/Zhang ZM/Lambrecht ====
# 1) Is there Tfh cluster? If so, heatmap showing signature genes of TFH. Especially IL21
# ASCL2, BCL6, CD4, CD200, CXCR5, IL4, IL21, IL6ST, MAF, PDCD1, SH2D1A, TOX2
# 2) Is there GC B cluster? If so, heatmap showing signature genes of GCB signature,
# AICDA, BATF, BACH2, BCL6, CD79A, CD79B, CD86, DOCK8, IRF4, IRF8, MYC

### plot data Kim.(w or w/o seurat) ====
## refine
load('kim_umi_lu.RData')
rm(meta_lu)
class(mat_lu) <- class(as.data.frame(mat_lu))
genes <- mat_lu$Index
mat_lu <- mat_lu[,-1]

annot_lu$index_new <- NA
for (num in levels(factor(annot_lu$number))){
  id <- which(annot_lu$number==num)
  annot_lu$index_new[id] <- paste0(num, "_", 1:length(id))
}
names(mat_lu) <- annot_lu$index_new
save(annot_lu, genes, mat_lu, file="kim_umi_lu.RData")


### ====



