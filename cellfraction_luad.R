#### Tfh/B cell proportion analysis 200510
setwd('~/Documents/Projects/cc/tfh_fraction/')
library(dplyr)

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



### END ###### 



