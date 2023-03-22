library( "DESeq2" )
library(ggplot2)
library(apeglm)
library(scales)
library(viridis)
library(devtools)
library(ggplot2)
library(PCAtools)
library(Cairo)
library(clusterProfiler)
library(gage)

## Tree B3
deseqFile1 <- "b3_counts.csv"
countData1 <- read.table(deseqFile1, header = T, sep = ";", row.names = 1)
head(countData1)
sampleNames1 <- colnames(countData1)
sampleCondition1 <- c("Control","Control","Control","Treatment", "Treatment","Treatment")
sampleCondition2 <- c("CH20","CH48","CH108","CH17", "CH45","CH105")

colData1<- data.frame(condition = sampleCondition1, tree=sampleCondition2)
row.names(colData1) = sampleNames1
treatments = c("Control","Treatment")
all(rownames(colData1) == colnames(countData1))
dds1 <- DESeqDataSetFromMatrix(countData = countData1,
                               colData = colData1,
                               design = ~ condition)
dds1
keep1 <- rowSums(counts(dds1)) >= 0
dds1 <- dds1[keep1,]

dds1 <- DESeq(dds1)
res1 <- results(dds1)
res1
res1 <- res1[!is.na(res1$padj),]


res1 <- results(dds1, name="condition_Treatment_vs_Control")
res1 <- results(dds1, contrast=c("condition", "Treatment","Control"))
resultsNames(dds1)
resLFC_1 <- lfcShrink(dds1, coef="condition_Treatment_vs_Control", type="apeglm")
head(resLFC_1)
resOrdered_1 <- resLFC_1[order(resLFC_1$pvalue),]
summary(resLFC_1)


resLFC_1$ID<-rownames(resLFC_1)

head(dds1)

sum(resLFC_1$padj < 0.1, na.rm=TRUE)
sum(resLFC_1$padj <= 0.05, na.rm=TRUE)

resOrdered <- res1[order(res1$pvalue),]
summary(res1)

sum(res1$padj < 0.1, na.rm=TRUE)
res05 <- results(dds1, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

plotMA(resLFC_1)
pvals<-as.data.frame(resOrdered_1[order(resOrdered_1$padj),])
sum(resLFC_1$padj <= 0.05 & pvals$log2FoldChange <= -2, na.rm=TRUE)
sum(resLFC_1$padj <= 0.05 & pvals$log2FoldChange >= 2, na.rm=TRUE)


vsd <- vst(dds1, blind=FALSE)

plotPCA(vsd, intgroup=c("condition", "tree"))
PCA_1<-plotPCA(vsd, intgroup="condition")+labs(title="B3")

library(EnhancedVolcano)

## Simple function for plotting a Volcano plot, returns a ggplot object
deseq.volcano <- function(res, datasetName) {
  return(EnhancedVolcano(res, x = 'log2FoldChange', y = 'padj',
                         #lab=rownames(res),
                         lab = NA,
                         #title = paste(datasetName, "Control vs EAB infested"),
                         title =NULL,
                         subtitle = NULL,
                         #subtitle = bquote(italic('FDR <= 0.05 and absolute FC >= 2')),
                         # Change text and icon sizes
                         labSize = 3, pointSize = 1.5, axisLabSize=10, titleLabSize=12,
                         subtitleLabSize=8, captionLabSize=10,
                         # Disable legend
                         legendPosition = "none",
                         # Set cutoffs
                         pCutoff = 0.05, FCcutoff = 2))
}

## Note: input data is the corrected DESeq2 output using the 'lfcShrink' function (see chapter 4)
p1<-deseq.volcano(res = resLFC_1, datasetName = "B3")

control_up_b3 <- pvals[pvals$padj < 0.05 & pvals$log2FoldChange < 0,]
treatment_up_b3 <- pvals[pvals$padj < 0.05 & pvals$log2FoldChange > 0,]
head(control_up_b3)
###########
## Tree B8
deseqFile1 <- "b8_counts.csv"
countData1 <- read.table(deseqFile1, header = T, sep = ";", row.names = 1)
head(countData1)
sampleNames1 <- colnames(countData1)
sampleCondition1 <- c("Control","Control","Control","Treatment", "Treatment","Treatment")
sampleCondition2 <- c("CH100","CH16","CH44","CH13", "CH41","CH97")

colData1<- data.frame(condition = sampleCondition1, tree=sampleCondition2)
row.names(colData1) = sampleNames1
treatments = c("Control","Treatment")
all(rownames(colData1) == colnames(countData1))
dds1 <- DESeqDataSetFromMatrix(countData = countData1,
                               colData = colData1,
                               design = ~ condition)
dds1
keep1 <- rowSums(counts(dds1)) >= 0
dds1 <- dds1[keep1,]

dds1 <- DESeq(dds1)
res1 <- results(dds1)

res1 <- res1[!is.na(res1$padj),]


res1 <- results(dds1, name="condition_Treatment_vs_Control")
res1 <- results(dds1, contrast=c("condition", "Treatment","Control"))
resultsNames(dds1)
resLFC_1 <- lfcShrink(dds1, coef="condition_Treatment_vs_Control", type="apeglm")
head(resLFC_1)
resOrdered_1 <- resLFC_1[order(resLFC_1$pvalue),]
summary(resLFC_1)


resLFC_1$ID<-rownames(resLFC_1)

head(dds1)

sum(resLFC_1$padj < 0.1, na.rm=TRUE)
sum(resLFC_1$padj <= 0.05 & resLFC_1$log2FoldChange <= -2, na.rm=TRUE)
sum(resLFC_1$padj <= 0.05 & resLFC_1$log2FoldChange >= 2, na.rm=TRUE)

vsd <- vst(dds1, blind=FALSE)

resOrdered <- res1[order(res1$pvalue),]
summary(res1)

sum(res1$padj < 0.1, na.rm=TRUE)
res05 <- results(dds1, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

#write.csv(as.data.frame(resOrdered_1[order(resOrdered_1$padj),] ), file="condition_control_vs_treatmentB8.csv")
plotMA(resLFC_1)
PCA_2<-plotPCA(vsd, intgroup="condition")+labs(title="B8")

pvals<-as.data.frame(resOrdered_1[order(resOrdered_1$padj),])
p2<-deseq.volcano(res = resLFC_1, datasetName = "B8")

control_up_b8 <- pvals[pvals$padj < 0.05 & pvals$log2FoldChange < 0,]
treatment_up_b8 <- pvals[pvals$padj < 0.05 & pvals$log2FoldChange > 0,]

###########
## Tree B9
deseqFile1 <- "b9_counts.csv"
countData1 <- read.table(deseqFile1, header = T, sep = ";", row.names = 1)
str(countData1)
sampleNames1 <- colnames(countData1)
sampleCondition1 <- c("Control","Control","Control","Treatment", "Treatment","Treatment")
sampleCondition2 <- c("CH4","CH32","CH88","CH1", "CH29","CH85")

colData1<- data.frame(condition = sampleCondition1, tree=sampleCondition2)
row.names(colData1) = sampleNames1
treatments = c("Control","Treatment")
all(rownames(colData1) == colnames(countData1))
dds1 <- DESeqDataSetFromMatrix(countData = countData1,
                               colData = colData1,
                               design = ~ condition)
dds1
keep1 <- rowSums(counts(dds1)) >= 0
dds1 <- dds1[keep1,]

dds1 <- DESeq(dds1)
res1 <- results(dds1)
res1
res1 <- res1[!is.na(res1$padj),]


res1 <- results(dds1, name="condition_Treatment_vs_Control")
res1 <- results(dds1, contrast=c("condition", "Treatment","Control"))
resultsNames(dds1)
resLFC_1 <- lfcShrink(dds1, coef="condition_Treatment_vs_Control", type="apeglm")
head(resLFC_1)
resOrdered_1 <- resLFC_1[order(resLFC_1$pvalue),]
summary(resLFC_1)


resLFC_1$ID<-rownames(resLFC_1)

head(dds1)

sum(resLFC_1$padj < 0.1, na.rm=TRUE)
sum(resLFC_1$padj < 0.05, na.rm=TRUE)
sum(resLFC_1$padj < 0.05 & resLFC_1$log2FoldChange <= -2, na.rm=TRUE)
sum(resLFC_1$padj < 0.05 & resLFC_1$log2FoldChange >= 2, na.rm=TRUE)

vsd <- vst(dds1, blind=FALSE)

resOrdered <- res1[order(res1$pvalue),]
summary(res1)

sum(res1$padj < 0.1, na.rm=TRUE)
res05 <- results(dds1, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

#write.csv(as.data.frame(resOrdered_1[order(resOrdered_1$padj),] ), file="condition_control_vs_treatmentB9.csv")
plotMA(resLFC_1)
pvals<-as.data.frame(resOrdered_1[order(resOrdered_1$padj),])
p3<-deseq.volcano(res = resLFC_1, datasetName = "B9")
PCA_3<-plotPCA(vsd, intgroup="condition")+labs(title="B9")


control_up_b9 <- pvals[pvals$padj < 0.05 & pvals$log2FoldChange < 0,]
treatment_up_b9 <- pvals[pvals$padj < 0.05 & pvals$log2FoldChange > 0,]
###########
## Tree B20
deseqFile1 <- "b20_counts.csv"
countData1 <- read.table(deseqFile1, header = T, sep = ",", row.names = 1)
#countData1[is.na(countData1)] <- 0

str(countData1)
sampleNames1 <- colnames(countData1)
sampleCondition1 <- c("Control","Control","Control","Treatment", "Treatment","Treatment")
sampleCondition2 <- c("CH8","CH36","CH64","CH5", "CH33","CH125")

colData1<- data.frame(condition = sampleCondition1, tree=sampleCondition2)
row.names(colData1) = sampleNames1
treatments = c("Control","Treatment")
all(rownames(colData1) == colnames(countData1))
dds1 <- DESeqDataSetFromMatrix(countData = countData1,
                               colData = colData1,
                               design = ~ condition)
dds1
keep1 <- rowSums(counts(dds1)) >= 0
dds1 <- dds1[keep1,]

dds1 <- DESeq(dds1)
res1 <- results(dds1)

res1 <- res1[!is.na(res1$padj),]


res1 <- results(dds1, name="condition_Treatment_vs_Control")
res1 <- results(dds1, contrast=c("condition", "Treatment","Control"))
resultsNames(dds1)
resLFC_1 <- lfcShrink(dds1, coef="condition_Treatment_vs_Control", type="apeglm")
head(resLFC_1)
resOrdered_1 <- resLFC_1[order(resLFC_1$pvalue),]
summary(resLFC_1)


resLFC_1$ID<-rownames(resLFC_1)

head(dds1)

sum(resLFC_1$padj < 0.1, na.rm=TRUE)
sum(resLFC_1$padj < 0.05, na.rm=TRUE)
sum(resLFC_1$padj <= 0.05 & resLFC_1$log2FoldChange >= 2, na.rm=TRUE)
sum(resLFC_1$padj <= 0.05 & resLFC_1$log2FoldChange <= -2, na.rm=TRUE)


vsd <- vst(dds1, blind=FALSE)

resOrdered <- res1[order(res1$pvalue),]
summary(res1)

sum(res1$padj < 0.1, na.rm=TRUE)
res05 <- results(dds1, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

#write.csv(as.data.frame(resOrdered_1[order(resOrdered_1$padj),] ), file="condition_control_vs_treatmentB20.csv")
plotMA(resLFC_1)
pvals<-as.data.frame(resOrdered_1[order(resOrdered_1$padj),])
p4<-deseq.volcano(res = resLFC_1, datasetName = "B20")
p4 +
  ggplot2::coord_cartesian(xlim=c(-10, 20), ylim=c(0, 60))



##########################
PCA_4<-plotPCA(vsd, intgroup="condition", returnData = FALSE)+labs(title="B20")

control_up_b20 <- pvals[pvals$padj < 0.05 & pvals$log2FoldChange < 0,]
treatment_up_b20 <- pvals[pvals$padj < 0.05 & pvals$log2FoldChange > 0,]
#write.csv(control_up_b20, "B20_up.csv")

######################

library(gridExtra)
library(grid)
grid.arrange(p1, p2,p3,p4,
             ncol=2,
             nrow=2,
             top = textGrob('EnhancedVolcano',
                            just = c('center'),
                            gp = gpar(fontsize = 32)))

p5<-grid.arrange(p1, p2,p3,p4,
                 ncol=2,
                 nrow=2)

png("Volcano_EAB.png", width = 10, height = 12, units = 'in', res = 300)
plot (p5)# Make plot
dev.off() 

svg("Volcano_EAB.svg", width = 7, height = 6)
plot (p5)# Make plot
dev.off() 
####################
# combined pca
grid.newpage()


PCA_5<-grid.arrange(PCA_1, PCA_2,PCA_3,PCA_4,
                    ncol=2,
                    nrow=2)

svg("PCA_suc_v_tol.svg", width = 7, height = 6)
plot (PCA_5)# Make plot
dev.off() 
###############
# combined venn diagram

# move to new plotting page
grid.newpage()

control_up_b3_genes <- row.names(control_up_b3)
control_up_b8_genes <- row.names(control_up_b8)
control_up_b9_genes <- row.names(control_up_b9)
control_up_b20_genes <- row.names(control_up_b20)
head(control_up_b3_genes)
treatment_up_b3_genes <- row.names(treatment_up_b3)
treatment_up_b8_genes <- row.names(treatment_up_b8)
treatment_up_b9_genes <- row.names(treatment_up_b9)
treatment_up_b20_genes <- row.names(treatment_up_b20)

library(gplots)

my_list <- list(control_up_b3_genes, control_up_b8_genes, control_up_b9_genes, control_up_b20_genes)
venn(my_list)
my_list2 <- list(
  B3 = control_up_b3_genes, 
  B8 = control_up_b8_genes, 
  B9 = control_up_b9_genes,
  B20 = control_up_b20_genes
)

my_list3 <- list(
  B3 = treatment_up_b3_genes, 
  B8 = treatment_up_b8_genes, 
  B9 = treatment_up_b9_genes,
  B20 = treatment_up_b20_genes
)
# save object from venn
my_venn <- venn(my_list2, show.plot = FALSE)
class(my_venn)

library(ggvenn)
q2<-ggvenn(
  my_list3, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

png("EAB_treatment_venn.png", width = 12, height = 8, units = 'in', res = 300)
plot (q2)# Make plot
dev.off() 

svg("EAB_treatment_venn.svg", width = 7, height = 6)
plot (q2)# Make plot
dev.off() 

library("ggVennDiagram")
q4<-ggVennDiagram(my_list3, label_alpha = 0) +
  scale_fill_gradient(low="blue",high = "red")

png("EAB_treatment_venn_heat.png", width = 12, height = 8, units = 'in', res = 300)
plot (q4)# Make plot
dev.off() 

svg("EAB_treatment_venn_heat.svg", width = 7, height = 6)
plot (q4)# Make plot
dev.off() 

sessionInfo()
###########
# Trees B3 v B8 

deseqFile1 <- "B3_v_B8.csv"
countData1 <- read.table(deseqFile1, header = T, sep = ",", row.names = 1)
head(countData1)
#countData1[is.na(countData1)] <- 0
sampleNames1 <- colnames(countData1)
sampleCondition1 <- c("Control","Control","Control","Treatment", "Treatment","Treatment", "Control", "Control", "Control", "Treatment", "Treatment", "Treatment")

#sampleCondition2 <- c("Suceptible","Suceptible","Suceptible", "Suceptible","Suceptible","Suceptible","Suceptible","Suceptible","Suceptible", "Suceptible","Suceptible","Suceptible", "Tolerant", "Tolerant", "Tolerant", "Tolerant", "Tolerant", "Tolerant", "Tolerant", "Tolerant", "Tolerant", "Tolerant", "Tolerant", "Tolerant")
sampleCondition2<- c("A","A","A","A","A","A","B","B","B","B","B","B")
colData1<- data.frame(condition=sampleCondition1, genotype=sampleCondition2)

colData1$condition<-as.factor(colData1$condition)
colData1$type<-as.factor(colData1$type)
colData1$genotype<-as.factor(colData1$genotype)

str(colData1)

row.names(colData1) = sampleNames1
#treatments = c("Control_Suceptible_1","Infested_Tolerant_1","Control_Suceptible_2","Infested_Tolerant_2")
all(rownames(colData1) == colnames(countData1))
dds1 <- DESeqDataSetFromMatrix(countData = countData1,
                               colData = colData1,
                               design = ~ condition)

colData(dds1)


dds <- DESeq(dds1)

resultsNames(dds)

res1 <- results(dds, name="condition_Treatment_vs_Control")
res1 <- results(dds, contrast=c("condition", "Treatment","Control"))
resultsNames(dds)
resLFC_1 <- lfcShrink(dds, coef="condition_Treatment_vs_Control", type="apeglm")
head(resLFC_1)
resOrdered_1 <- resLFC_1[order(resLFC_1$pvalue),]
summary(resLFC_1)

res1 <- res1[!is.na(res1$padj),]

resLFC_1$ID<-rownames(resLFC_1)

head(dds)

sum(resLFC_1$padj < 0.1, na.rm=TRUE)
sum(resLFC_1$padj <= 0.05, na.rm=TRUE)

resOrdered <- res1[order(res1$pvalue),]
summary(res1)

sum(res1$padj < 0.1, na.rm=TRUE)
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

vsd <- vst(dds, blind=FALSE)

resOrdered <- res1[order(res1$pvalue),]
summary(res1)

sum(res1$padj < 0.1, na.rm=TRUE)
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

#write.csv(as.data.frame(resOrdered_1[order(resOrdered_1$padj),] ), file="condition_b3_b8_vs_b9_b20_treatment_only.csv")
plotMA(resLFC_1, ylim=c(-8,8), cex=.8)
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
plotDispEsts(dds)

pvals<-as.data.frame(resOrdered_1[order(resOrdered_1$padj),])
p6<-deseq.volcano(res = resLFC_1, datasetName = "B3 v B8")
#p4 +
#  ggplot2::coord_cartesian(xlim=c(-10, 20), ylim=c(0, 60))

png("Volcano_EAB_B3_v_B8.png", width = 10, height = 12, units = 'in', res = 300)
plot (p6)# Make plot
dev.off() 

svg("Volcano_EAB_B3_v_B8.svg", width = 7, height = 6)
plot (p6)# Make plot
dev.off() 


PCA_5<-plotPCA(vsd, intgroup="condition")+labs(title="B3 v B8")
png("PCA_B3_v_B8.png", width = 10, height = 12, units = 'in', res = 300)
plot (PCA_5)# Make plot
dev.off() 
svg("PCA_B3_v_B8.svg", width = 7, height = 6)
plot (PCA_5)# Make plot
dev.off() 

control_up_suceptible <- pvals[pvals$padj < 0.05 & pvals$log2FoldChange < 0,]
treatment_up_tolerant <- pvals[pvals$padj < 0.05 & pvals$log2FoldChange > 0,]

control_up_suceptible <- row.names(control_up_suceptible)
head(control_up_suceptible)
treatment_up_tolerant <- row.names(treatment_up_tolerant)

library(gplots)

my_list <- list(control_up_suceptible,treatment_up_tolerant)
venn(my_list)
my_list2 <- list(
  control = control_up_suceptible, 
  infested = treatment_up_tolerant
)

my_venn <- venn(my_list2, show.plot = FALSE)
class(my_venn)

library(ggvenn)
q2<-ggvenn(
  my_list2, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

png("EAB_B3_v_B8_venn.png", width = 12, height = 8, units = 'in', res = 300)
plot (q2)# Make plot
dev.off() 

svg("EAB_B3_v_B8_venn.svg", width = 7, height = 6)
plot (q2)# Make plot
dev.off() 

###########
# Trees B9 v B20 

deseqFile1 <- "B9_v_B20.csv"
countData1 <- read.table(deseqFile1, header = T, sep = ",", row.names = 1)
head(countData1)
#countData1[is.na(countData1)] <- 0
sampleNames1 <- colnames(countData1)
sampleCondition1 <- c("Control","Control","Control","Treatment", "Treatment","Treatment", "Control", "Control", "Control", "Treatment", "Treatment", "Treatment")

#sampleCondition2 <- c("Suceptible","Suceptible","Suceptible", "Suceptible","Suceptible","Suceptible","Suceptible","Suceptible","Suceptible", "Suceptible","Suceptible","Suceptible", "Tolerant", "Tolerant", "Tolerant", "Tolerant", "Tolerant", "Tolerant", "Tolerant", "Tolerant", "Tolerant", "Tolerant", "Tolerant", "Tolerant")
#sampleCondition2<- c("C","C","C","A","A","A","B","B","B","B","B","B")
colData1<- data.frame(condition=sampleCondition1)

colData1$condition<-as.factor(colData1$condition)

str(colData1)

row.names(colData1) = sampleNames1
#treatments = c("Control_Suceptible_1","Infested_Tolerant_1","Control_Suceptible_2","Infested_Tolerant_2")
all(rownames(colData1) == colnames(countData1))
dds1 <- DESeqDataSetFromMatrix(countData = countData1,
                               colData = colData1,
                               design = ~ condition)

colData(dds1)


dds <- DESeq(dds1)

resultsNames(dds)

res1 <- results(dds, name="condition_Treatment_vs_Control")
res1 <- results(dds, contrast=c("condition", "Treatment","Control"))
resultsNames(dds)
resLFC_1 <- lfcShrink(dds, coef="condition_Treatment_vs_Control", type="apeglm")
head(resLFC_1)
resOrdered_1 <- resLFC_1[order(resLFC_1$pvalue),]
summary(resLFC_1)

res1 <- res1[!is.na(res1$padj),]

resLFC_1$ID<-rownames(resLFC_1)

head(dds)

sum(resLFC_1$padj < 0.1, na.rm=TRUE)
sum(resLFC_1$padj <= 0.05, na.rm=TRUE)

resOrdered <- res1[order(res1$pvalue),]
summary(res1)

sum(res1$padj < 0.1, na.rm=TRUE)
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

vsd <- vst(dds, blind=FALSE)

resOrdered <- res1[order(res1$pvalue),]
summary(res1)

sum(res1$padj < 0.1, na.rm=TRUE)
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

#write.csv(as.data.frame(resOrdered_1[order(resOrdered_1$padj),] ), file="condition_b3_b8_vs_b9_b20_treatment_only.csv")
plotMA(resLFC_1, ylim=c(-8,8), cex=.8)
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
plotDispEsts(dds)

pvals<-as.data.frame(resOrdered_1[order(resOrdered_1$padj),])
p7<-deseq.volcano(res = resLFC_1, datasetName = "B9 v B20")
#p4 +
#  ggplot2::coord_cartesian(xlim=c(-10, 20), ylim=c(0, 60))

png("Volcano_EAB_B9_v_B20.png", width = 10, height = 12, units = 'in', res = 300)
plot (p7)# Make plot
dev.off() 

svg("Volcano_EAB_B9_v_B20.svg", width = 7, height = 6)
plot (p7)# Make plot
dev.off() 


PCA_6<-plotPCA(vsd, intgroup="condition")+labs(title="B9 v B20")
png("PCA_B9_v_B20.png", width = 10, height = 12, units = 'in', res = 300)
plot (PCA_6)# Make plot
dev.off() 
svg("PCA_B9_v_B20.svg", width = 7, height = 6)
plot (PCA_6)# Make plot
dev.off() 

control_up_suceptible <- pvals[pvals$padj < 0.05 & pvals$log2FoldChange < 0,]
treatment_up_tolerant <- pvals[pvals$padj < 0.05 & pvals$log2FoldChange > 0,]

control_up_suceptible <- row.names(control_up_suceptible)
head(control_up_suceptible)
treatment_up_tolerant <- row.names(treatment_up_tolerant)

library(gplots)

my_list <- list(control_up_suceptible,treatment_up_tolerant)
venn(my_list)
my_list2 <- list(
  control = control_up_suceptible, 
  infested = treatment_up_tolerant
)

my_venn <- venn(my_list2, show.plot = FALSE)
class(my_venn)

library(ggvenn)
q3<-ggvenn(
  my_list2, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

png("EAB_B9_v_B20_venn.png", width = 12, height = 8, units = 'in', res = 300)
plot (q3)# Make plot
dev.off() 

svg("EAB_B9_v_B20_venn.svg", width = 7, height = 6)
plot (q3)# Make plot
dev.off() 
########################
# Trees all clone 

deseqFile1 <- "B3_B8_v_B9_B20_treatment_clean.csv"
countData1 <- read.table(deseqFile1, header = T, sep = ",", row.names = 1)
head(countData1)
#countData1[is.na(countData1)] <- 0
sampleNames1 <- colnames(countData1)
sampleCondition1 <- c("Control","Control","Control","Treatment", "Treatment","Treatment", "Control", "Control", "Control", "Treatment", "Treatment", "Treatment","Control","Control","Control","Treatment", "Treatment","Treatment", "Control", "Control", "Control", "Treatment", "Treatment", "Treatment")

sampleCondition2 <- c("Control_suc_1","Control_suc_1","Control_suc_1", "Infested_suc_1","Infested_suc_1","Infested_suc_1","Control_suc_2","Control_suc_2","Control_suc_2", "Infested_suc_2","Infested_suc_2","Infested_suc_2", "Control_tol_1", "Control_tol_1", "Control_tol_1", "Infested_tol_1", "Infested_tol_1", "Infested_tol_1", "Control_tol_2", "Control_tol_2", "Control_tol_2", "Infested_tol_2", "Infested_tol_2", "Infested_tol_2")
colData1<- data.frame(condition=sampleCondition1, four_groups=sampleCondition2)

colData1$condition<-as.factor(colData1$condition)
colData1$four_groups<-as.factor(colData1$four_groups)

str(colData1)

row.names(colData1) = sampleNames1
#treatments = c("Control_Suceptible_1","Infested_Tolerant_1","Control_Suceptible_2","Infested_Tolerant_2")
all(rownames(colData1) == colnames(countData1))
dds1 <- DESeqDataSetFromMatrix(countData = countData1,
                               colData = colData1,
                               design = ~ condition)

colData(dds1)


dds <- DESeq(dds1)

resultsNames(dds)

res1 <- results(dds, name="condition_Treatment_vs_Control")
res1 <- results(dds, contrast=c("condition", "Treatment","Control"))
resultsNames(dds)
resLFC_1 <- lfcShrink(dds, coef="condition_Treatment_vs_Control", type="apeglm")
head(resLFC_1)
resOrdered_1 <- resLFC_1[order(resLFC_1$pvalue),]
summary(resLFC_1)

res1 <- res1[!is.na(res1$padj),]

resLFC_1$ID<-rownames(resLFC_1)

head(dds)

sum(resLFC_1$padj < 0.1, na.rm=TRUE)
sum(resLFC_1$padj <= 0.05, na.rm=TRUE)

resOrdered <- res1[order(res1$pvalue),]
summary(res1)

sum(res1$padj < 0.1, na.rm=TRUE)
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

vsd <- vst(dds, blind=FALSE)

resOrdered <- res1[order(res1$pvalue),]
summary(res1)

sum(res1$padj < 0.1, na.rm=TRUE)
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

#write.csv(as.data.frame(resOrdered_1[order(resOrdered_1$padj),] ), file="condition_b3_b8_vs_b9_b20_treatment_only.csv")
plotMA(resLFC_1, ylim=c(-8,8), cex=.8)
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
plotDispEsts(dds)

pvals<-as.data.frame(resOrdered_1[order(resOrdered_1$padj),])
p8<-deseq.volcano(res = resLFC_1, datasetName = "all_samples")
#p4 +
#  ggplot2::coord_cartesian(xlim=c(-10, 20), ylim=c(0, 60))

png("Volcano_control_v_treatment_all_samples.png", width = 10, height = 12, units = 'in', res = 300)
plot (p8)# Make plot
dev.off() 

svg("Volcano_control_v_treatment_all_samples.svg", width = 7, height = 6)
plot (p8)# Make plot
dev.off() 


PCA_7<-plotPCA(vsd, intgroup="four_groups")+labs(title="all_samples")
png("PCA_all_2.png", width = 10, height = 12, units = 'in', res = 300)
plot (PCA_7)# Make plot
dev.off() 
svg("PCA_all.svg", width = 7, height = 6)
plot (PCA_7)# Make plot
dev.off() 

################
# enrichment
library(KEGGREST)
library(gage)

##
# get data from deseq2

pvals<-as.data.frame(resOrdered_1[order(resOrdered_1$padj),])
pvals$gene <- rownames(pvals)

# these annotation files are on github
olea<-read.csv("B3_blast_out4.csv", sep = "\t", header = FALSE)
olea_2<-read.csv("proteins_10724_352035.csv", sep = ",", header = TRUE)
olea$gi<-olea_2$GeneID[match(olea$V2, olea_2$Protein.product)]
olea3<-olea[, c(1,13)]
olea3$V1 = substr(olea3$V1, 1, nchar(olea3$V1)-2)

pvals$gi<-olea3$gi[match(pvals$gene, olea3$V1)]
control<- pvals[pvals$padj <= 0.05,]

pvals2<-control[, c(2,7)]
pvals3<-na.omit(pvals2)

#write.csv(pvals3, "B3_with_gi.csv")
pvals4<-read.csv("B3_with_gi.csv", header = TRUE)

str(pvals4)

pvals5<-pvals4[!duplicated(pvals4$gi), ]
rownames(pvals5) <- pvals5$gi

pvals6<-pvals5[, c(2,3)]

ref.idx=1
samp.idx=2
ms=apply(pvals6, 1, mean) 
exprs2=pvals6-ms

kg.nsy=kegg.gsets("nsy")
kg.nsy.eg=kegg.gsets("nsy", id.type = "kegg")

kg.cic=kegg.gsets("cic")
kg.cic.eg=kegg.gsets("cic", id.type = "kegg")

kg.egr=kegg.gsets("egr")
kg.egr.eg=kegg.gsets("egr", id.type = "kegg")

kg.nnu=kegg.gsets("nnu")
kg.nnu.eg=kegg.gsets("nnu", id.type = "kegg")

kg.pda=kegg.gsets("pda")
kg.pda.eg=kegg.gsets("pda", id.type = "kegg")

kg.rcu=kegg.gsets("rcu")
kg.rcu.eg=kegg.gsets("rcu", id.type = "kegg")

kg.pvu=kegg.gsets("pvu")
kg.pvu.eg=kegg.gsets("pvu", id.type = "kegg")

kg.mnt=kegg.gsets("mnt")
kg.mnt.eg=kegg.gsets("mnt", id.type = "kegg")

kg.cmo=kegg.gsets("cmo")
kg.cmo.eg=kegg.gsets("cmo", id.type = "kegg")

kg.cam=kegg.gsets("cam")
kg.cam.eg=kegg.gsets("cam", id.type = "kegg")

kg.bvg=kegg.gsets("bvg")
kg.bvg.eg=kegg.gsets("bvg", id.type = "kegg")

kg.crb=kegg.gsets("crb")
kg.crb.eg=kegg.gsets("crb", id.type = "kegg")

kg.cam=kegg.gsets("cam")
kg.cam.eg=kegg.gsets("cam", id.type = "kegg")

kg.nnu=kegg.gsets("nnu")
kg.nnu.eg=kegg.gsets("nnu", id.type = "kegg")

kg.oeu=kegg.gsets("oeu")
kg.oeu.eg=kegg.gsets("oeu", id.type = "kegg")

kg.ath=kegg.gsets("ath")
kg.ath.eg=kegg.gsets("ath", id.type = "kegg")

combined_kg_eg_sets=c(kg.nsy.eg$kg.sets, kg.nsy.eg$sigmet.idx, kg.nsy.eg$sig.idx, kg.nsy.eg$met.idx, kg.nsy.eg$dise.idx, kg.cic.eg$kg.sets, kg.cic.eg$sigmet.idx, kg.cic.eg$sig.idx, kg.cic.eg$met.idx, kg.cic.eg$dise.idx, kg.egr.eg$kg.sets, kg.egr.eg$sigmet.idx, kg.egr.eg$sig.idx, kg.egr.eg$met.idx, kg.egr.eg$dise.idx,kg.nnu.eg$kg.sets, kg.nnu.eg$sigmet.idx, kg.nnu.eg$sig.idx, kg.nnu.eg$met.idx, kg.nnu.eg$dise.idx,kg.pda.eg$kg.sets, kg.pda.eg$sigmet.idx, kg.pda.eg$sig.idx, kg.pda.eg$met.idx, kg.pda.eg$dise.idx,kg.rcu.eg$kg.sets, kg.rcu.eg$sigmet.idx, kg.rcu.eg$sig.idx, kg.rcu.eg$met.idx, kg.rcu.eg$dise.idx,kg.pvu.eg$kg.sets, kg.pvu.eg$sigmet.idx, kg.pvu.eg$sig.idx, kg.pvu.eg$met.idx, kg.pvu.eg$dise.idx,kg.mnt.eg$kg.sets, kg.mnt.eg$sigmet.idx, kg.mnt.eg$sig.idx, kg.mnt.eg$met.idx, kg.mnt.eg$dise.idx,kg.cmo.eg$kg.sets, kg.cmo.eg$sigmet.idx, kg.cmo.eg$sig.idx, kg.cmo.eg$met.idx, kg.cmo.eg$dise.idx,kg.cam.eg$kg.sets, kg.cam.eg$sigmet.idx, kg.cam.eg$sig.idx, kg.cam.eg$met.idx, kg.cam.eg$dise.idx,kg.bvg.eg$kg.sets, kg.bvg.eg$sigmet.idx, kg.bvg.eg$sig.idx, kg.bvg.eg$met.idx, kg.bvg.eg$dise.idx,kg.crb.eg$kg.sets, kg.crb.eg$sigmet.idx, kg.crb.eg$sig.idx, kg.crb.eg$met.idx, kg.crb.eg$dise.idx, kg.oeu.eg$kg.sets, kg.oeu.eg$sigmet.idx, kg.oeu.eg$sig.idx, kg.oeu.eg$met.idx, kg.oeu.eg$dise.idx, kg.ath.eg$kg.sets, kg.ath.eg$sigmet.idx, kg.ath.eg$sig.idx, kg.ath.eg$met.idx, kg.ath.eg$dise.idx)
oeu_only=c(kg.oeu.eg$kg.sets, kg.oeu.eg$sigmet.idx, kg.oeu.eg$sig.idx, kg.oeu.eg$met.idx, kg.oeu.eg$dise.idx)
keggresA = gage(exprs2, gsets=combined_kg_eg_sets,ref = samp.idx,  samp = ref.idx, same.dir   = F, compare = 'paired', set.size=c(50, 1000))
lapply(keggresA, head)

keggres2<-as.data.frame(keggresA)
head(keggres2)
keggres2$names <- rownames(keggres2)
keggres3<-transform(keggres2, KeggID = substr(names, 1,8), Annotations = substr(names, 1,40))
keggres4<-na.omit(keggres3)
Kdf5<-subset(keggres4, greater.p.val < 0.5) 

#write.csv(Kdf5, "Gene family enrichment B3 control v B3 ADB")

(q1 <- ggplot(Kdf5, aes(greater.stat.mean, Annotations))+
    geom_point(aes(colour=greater.q.val, size=greater.set.size)) +
    scale_size(name = "Number of genes")+
    ggtitle("B3 control vs EAB infested") +
    scale_color_gradientn(colours=rainbow(4), limits=c(0, 1), name = "q-value") +
    #geom_vline(xintercept=0, size=0.5, colour="gray50") +
    xlab("Magnitude of gene-set level changes")+
    labs(fill = "Number of genes")+
    #facet_wrap(~ assay_name, ncol=2)+
    theme(panel.background=element_rect(fill="NA", colour="black", size = 2, linetype='solid'),
          panel.grid.major=element_line(size=0.05,linetype='solid', colour="gray90"), 
          panel.grid.minor=element_line(size=0.05,linetype='solid', colour="gray90"),
          axis.text.y = element_text(size = 6),
          strip.text.x = element_text(size = 7),
          axis.title.y=element_blank()) +
    expand_limits(x=c(0, 5)) +
    scale_x_continuous(breaks=c(0,1,2,3,4)))
#legend.position = "none"))

png("Enrichment_genes_in_B3.png", res = 300, width = 2500, height = 1600)
plot (q1)# Make plot
dev.off() 

svg("Enrichment_genes_in_B3.svg",width = 7, height = 6)
plot (q1)# Make plot
dev.off() 
#################
# B8
#### from deseq2 in enrichment_EAB.R script
pvals<-as.data.frame(resOrdered_1[order(resOrdered_1$padj),])
pvals$gene <- rownames(pvals)
###
olea<-read.csv("B3_blast_out4.csv", sep = "\t", header = FALSE)
olea_2<-read.csv("proteins_10724_352035.csv", sep = ",", header = TRUE)
olea$gi<-olea_2$GeneID[match(olea$V2, olea_2$Protein.product)]
olea3<-olea[, c(1,13)]
olea3$V1 = substr(olea3$V1, 1, nchar(olea3$V1)-2)

pvals$gi<-olea3$gi[match(pvals$gene, olea3$V1)]
control<- pvals[pvals$padj <= 0.05,]

pvals2<-control[, c(2,7)]
pvals3<-na.omit(pvals2)

write.csv(pvals3, "B8_with_gi.csv")
pvals4<-read.csv("B8_with_gi.csv", header = TRUE)
str(pvals4)

pvals5<-pvals4[!duplicated(pvals4$gi), ]
rownames(pvals5) <- pvals5$gi
pvals6<-pvals5[, c(2,3)]

ref.idx=1
samp.idx=2
ms=apply(pvals6, 1, mean) 
exprs2=pvals6-ms

keggresA = gage(exprs2, gsets=combined_kg_eg_sets,ref = ref.idx,  samp = samp.idx, same.dir   = F, compare = 'paired', set.size=c(50, 1000))
lapply(keggresA, head)

keggres2<-as.data.frame(keggresA)
head(keggres2)
keggres2$names <- rownames(keggres2)
keggres3<-transform(keggres2, KeggID = substr(names, 1,8), Annotations = substr(names, 1,40))
keggres4<-na.omit(keggres3)
Kdf5<-subset(keggres4, greater.p.val < 0.5) 

#write.csv(Kdf5, "Gene family enrichment B8 control v B3 EAB")

(q2 <- ggplot(Kdf5, aes(greater.stat.mean, Annotations))+
    geom_point(aes(colour=greater.q.val, size=greater.set.size)) +
    scale_size(name = "Number of genes")+
    ggtitle("B8 control v EAB infested") +
    scale_color_gradientn(colours=rainbow(4), limits=c(0, 1), name = "q-value") +
    #geom_vline(xintercept=0, size=0.5, colour="gray50") +
    xlab("Magnitude of gene-set level changes")+
    labs(fill = "Number of genes")+
    #facet_wrap(~ assay_name, ncol=2)+
    theme(panel.background=element_rect(fill="NA", colour="black", size = 2, linetype='solid'),
          panel.grid.major=element_line(size=0.05,linetype='solid', colour="gray90"), 
          panel.grid.minor=element_line(size=0.05,linetype='solid', colour="gray90"),
          axis.text.y = element_text(size = 6),
          strip.text.x = element_text(size = 7),
          axis.title.y=element_blank()) +
    expand_limits(x=c(0, 5)) +
    scale_x_continuous(breaks=c(0,1,2,3,4,5)))
#legend.position = "none"))

png("Enrichment_genes_in_B8.png", res = 300, width = 2500, height = 1600)
plot (q2)# Make plot
dev.off() 

svg("Enrichment_genes_in_B8.svg",width = 7, height = 6)
plot (q2)# Make plot
dev.off() 
##############
## B9

pvals<-as.data.frame(resOrdered_1[order(resOrdered_1$padj),])
pvals$gene <- rownames(pvals)
###
olea<-read.csv("B3_blast_out4.csv", sep = "\t", header = FALSE)
olea_2<-read.csv("proteins_10724_352035.csv", sep = ",", header = TRUE)
olea$gi<-olea_2$GeneID[match(olea$V2, olea_2$Protein.product)]
olea3<-olea[, c(1,13)]
olea3$V1 = substr(olea3$V1, 1, nchar(olea3$V1)-2)

pvals$gi<-olea3$gi[match(pvals$gene, olea3$V1)]
control<- pvals[pvals$padj <= 0.05,]

pvals2<-control[, c(2,7)]
pvals3<-na.omit(pvals2)

#write.csv(pvals3, "B9_with_gi.csv")
pvals4<-read.csv("B9_with_gi.csv", header = TRUE)
str(pvals4)

pvals5<-pvals4[!duplicated(pvals4$gi), ]
rownames(pvals5) <- pvals5$gi
pvals6<-pvals5[, c(2,3)]

ref.idx=1
samp.idx=2
ms=apply(pvals6, 1, mean) 
exprs2=pvals6-ms

keggresA = gage(exprs2, gsets=combined_kg_eg_sets,ref = ref.idx,  samp = samp.idx, same.dir   = F, compare = 'paired', set.size=c(50, 1000))
lapply(keggresA, head)

keggres2<-as.data.frame(keggresA)
head(keggres2)
keggres2$names <- rownames(keggres2)
keggres3<-transform(keggres2, KeggID = substr(names, 1,8), Annotations = substr(names, 1,40))
keggres4<-na.omit(keggres3)
Kdf5<-subset(keggres4, greater.p.val < 0.5) 

#write.csv(Kdf5, "Gene family enrichment B9 control v B9 EAB")

(q3 <- ggplot(Kdf5, aes(greater.stat.mean, Annotations))+
    geom_point(aes(colour=greater.q.val, size=greater.set.size)) +
    scale_size(name = "Number of genes")+
    ggtitle("B9 control v EAB infested") +
    scale_color_gradientn(colours=rainbow(4), limits=c(0, 1), name = "q-value") +
    #geom_vline(xintercept=0, size=0.5, colour="gray50") +
    xlab("Magnitude of gene-set level changes")+
    labs(fill = "Number of genes")+
    #facet_wrap(~ assay_name, ncol=2)+
    theme(panel.background=element_rect(fill="NA", colour="black", size = 2, linetype='solid'),
          panel.grid.major=element_line(size=0.05,linetype='solid', colour="gray90"), 
          panel.grid.minor=element_line(size=0.05,linetype='solid', colour="gray90"),
          axis.text.y = element_text(size = 6),
          strip.text.x = element_text(size = 7),
          axis.title.y=element_blank()) +
    expand_limits(x=c(0, 5)) +
    scale_x_continuous(breaks=c(0,1,2,3,4,5)))
#legend.position = "none"))

png("Enrichment_genes_in_B9.png", res = 300, width = 2500, height = 1600)
plot (q3)# Make plot
dev.off() 

svg("Enrichment_genes_in_B9.svg",width = 7, height = 6)
plot (q3)# Make plot
dev.off() 

####################
## B20

pvals<-as.data.frame(resOrdered_1[order(resOrdered_1$padj),])
pvals$gene <- rownames(pvals)
###
olea<-read.csv("../B3_blast_out4.csv", sep = "\t", header = FALSE)
olea_2<-read.csv("../proteins_10724_352035.csv", sep = ",", header = TRUE)
olea$gi<-olea_2$GeneID[match(olea$V2, olea_2$Protein.product)]
olea3<-olea[, c(1,13)]
olea3$V1 = substr(olea3$V1, 1, nchar(olea3$V1)-2)

pvals$gi<-olea3$gi[match(pvals$gene, olea3$V1)]
control<- pvals[pvals$padj <= 0.05,]

pvals2<-control[, c(2,7)]
pvals3<-na.omit(pvals2)

write.csv(pvals3, "suc_tole_all_with_gi.csv")
pvals4<-read.csv("suc_tole_all_with_gi.csv", header = TRUE)
str(pvals4)

pvals5<-pvals4[!duplicated(pvals4$gi), ]
rownames(pvals5) <- pvals5$gi
pvals6<-pvals5[, c(2,3)]

ref.idx=1
samp.idx=2
ms=apply(pvals6, 1, mean) 
exprs2=pvals6-ms

keggresA = gage(pvals6, gsets=combined_kg_eg_sets,ref = ref.idx,  samp = samp.idx, same.dir   = F, compare = 'paired', set.size=c(50, 1000))
lapply(keggresA, head)

keggres2<-as.data.frame(keggresA)
head(keggres2)
keggres2$names <- rownames(keggres2)
keggres3<-transform(keggres2, names = substr(names, 1,8), Annotations = substr(names, 1,40))
keggres4<-na.omit(keggres3)
Kdf5<-subset(keggres4, greater.p.val < 0.5) 

#write.csv(Kdf5, "Gene family enrichment B20 control v B20 EAB")

(q4 <- ggplot(Kdf5, aes(greater.stat.mean, Annotations))+
    geom_point(aes(colour=greater.q.val, size=greater.set.size)) +
    scale_size(name = "Number of genes")+
    ggtitle("B20 control v EAB infested") +
    scale_color_gradientn(colours=rainbow(4), limits=c(0, 1), name = "q-value") +
    #geom_vline(xintercept=0, size=0.5, colour="gray50") +
    xlab("Magnitude of gene-set level changes")+
    labs(fill = "Number of genes")+
    #facet_wrap(~ assay_name, ncol=2)+
    theme(panel.background=element_rect(fill="NA", colour="black", size = 2, linetype='solid'),
          panel.grid.major=element_line(size=0.05,linetype='solid', colour="gray90"), 
          panel.grid.minor=element_line(size=0.05,linetype='solid', colour="gray90"),
          axis.text.y = element_text(size = 6),
          strip.text.x = element_text(size = 7),
          axis.title.y=element_blank()) +
    expand_limits(x=c(0, 5)) +
    scale_x_continuous(breaks=c(0,1,2,3,4,5)))
#legend.position = "none"))

png("Enrichment_genes_in_B20.png", res = 300, width = 2500, height = 1600)
plot (q4)# Make plot
dev.off() 

svg("Enrichment_genes_in_B20.svg",width = 7, height = 6)
plot (q4)# Make plot
dev.off() 
##############
library(gridExtra)
library(grid)
q5<-grid.arrange(q1, q2,q3,q4,
                 ncol=2,
                 nrow=2)
svg("Enrichment_EAB.svg", width = 7, height = 6)
plot (q5)# Make plot
dev.off() 

##############
## B3 vs B8

pvals<-as.data.frame(resOrdered_1[order(resOrdered_1$padj),])
pvals$gene <- rownames(pvals)

olea<-read.csv("../B3_blast_out4.csv", sep = "\t", header = FALSE)
olea_2<-read.csv("../proteins_10724_352035.csv", sep = ",", header = TRUE)
olea$gi<-olea_2$GeneID[match(olea$V2, olea_2$Protein.product)]
olea3<-olea[, c(1,13)]
olea3$V1 = substr(olea3$V1, 1, nchar(olea3$V1)-2)

pvals$gi<-olea3$gi[match(pvals$gene, olea3$V1)]
control<- pvals[pvals$padj <= 0.05,]

pvals2<-control[, c(2,7)]
pvals3<-na.omit(pvals2)

write.csv(pvals3, "B3_v_B8_gi.csv")
pvals4<-read.csv("B3_v_B8_gi.csv", header = TRUE)
str(pvals4)

pvals5<-pvals4[!duplicated(pvals4$gi), ]
rownames(pvals5) <- pvals5$gi
pvals6<-pvals5[, c(2,3)]

str(pvals6)
ref.idx=1
samp.idx=2
ms=apply(pvals6, 1, mean) 
exprs2=pvals6-ms
keggresA = gage(exprs2, gsets=combined_kg_eg_sets,ref = ref.idx,  samp = samp.idx, same.dir   = F, compare = 'paired', set.size=c(50, 1000))
lapply(keggresA, head)

keggres2<-as.data.frame(keggresA)
head(keggres2)
keggres2$names <- rownames(keggres2)
keggres3<-transform(keggres2, KeggID = substr(names, 1,8), Annotations = substr(names, 1,40))
keggres4<-na.omit(keggres3)
Kdf5<-subset(keggres4, greater.p.val < 0.5) 

#write.csv(Kdf5, "Gene family enrichment B3 control v B3 ADB")

(q5 <- ggplot(Kdf5, aes(greater.stat.mean, Annotations))+
    geom_point(aes(colour=greater.q.val, size=greater.set.size)) +
    scale_size(name = "Number of genes")+
    ggtitle("B3,B8 control v B3,B8 infested") +
    scale_color_gradientn(colours=rainbow(4), limits=c(0, 1), name = "q-value") +
    #geom_vline(xintercept=0, size=0.5, colour="gray50") +
    xlab("Magnitude of gene-set level changes")+
    labs(fill = "Number of genes")+
    #facet_wrap(~ assay_name, ncol=2)+
    theme(panel.background=element_rect(fill="NA", colour="black", size = 2, linetype='solid'),
          panel.grid.major=element_line(size=0.05,linetype='solid', colour="gray90"), 
          panel.grid.minor=element_line(size=0.05,linetype='solid', colour="gray90"),
          axis.text.y = element_text(size = 6),
          strip.text.x = element_text(size = 7),
          axis.title.y=element_blank()) +
    expand_limits(x=c(0, 5)) +
    scale_x_continuous(breaks=c(0,1,2,3,4,5)))
#legend.position = "none"))

png("Enrichment_genes_in_EAB_treatemnt only.png", res = 300, width = 2500, height = 1600)
plot (q5)# Make plot
dev.off() 

svg("Enrichment_B3VB8.svg", width = 7, height = 6)
plot (q5)# Make plot
dev.off() 

##############
## B9 vs B20

pvals<-as.data.frame(resOrdered_1[order(resOrdered_1$padj),])
pvals$gene <- rownames(pvals)

olea<-read.csv("../B3_blast_out4.csv", sep = "\t", header = FALSE)
olea_2<-read.csv("../proteins_10724_352035.csv", sep = ",", header = TRUE)
olea$gi<-olea_2$GeneID[match(olea$V2, olea_2$Protein.product)]
olea3<-olea[, c(1,13)]
olea3$V1 = substr(olea3$V1, 1, nchar(olea3$V1)-2)

pvals$gi<-olea3$gi[match(pvals$gene, olea3$V1)]
control<- pvals[pvals$padj <= 0.05,]

pvals2<-control[, c(2,7)]
pvals3<-na.omit(pvals2)

write.csv(pvals3, "B9_v_B20_gi.csv")
pvals4<-read.csv("B9_v_B20_gi.csv", header = TRUE)
str(pvals4)

pvals5<-pvals4[!duplicated(pvals4$gi), ]
rownames(pvals5) <- pvals5$gi
pvals6<-pvals5[, c(2,3)]

str(pvals6)
ref.idx=1
samp.idx=2
ms=apply(pvals6, 1, mean) 
exprs2=pvals6-ms
keggresA = gage(exprs2, gsets=combined_kg_eg_sets,ref = ref.idx,  samp = samp.idx, same.dir   = F, compare = 'paired', set.size=c(50, 1000))
lapply(keggresA, head)

keggres2<-as.data.frame(keggresA)
head(keggres2)
keggres2$names <- rownames(keggres2)
keggres3<-transform(keggres2, KeggID = substr(names, 1,8), Annotations = substr(names, 1,40))
keggres4<-na.omit(keggres3)
Kdf5<-subset(keggres4, greater.p.val < 0.5) 

#write.csv(Kdf5, "Gene family enrichment B3 control v B3 ADB")

(q6 <- ggplot(Kdf5, aes(greater.stat.mean, Annotations))+
    geom_point(aes(colour=greater.q.val, size=greater.set.size)) +
    scale_size(name = "Number of genes")+
    ggtitle("B9,B20 control v B9,B20 infested") +
    scale_color_gradientn(colours=rainbow(4), limits=c(0, 1), name = "q-value") +
    #geom_vline(xintercept=0, size=0.5, colour="gray50") +
    xlab("Magnitude of gene-set level changes")+
    labs(fill = "Number of genes")+
    #facet_wrap(~ assay_name, ncol=2)+
    theme(panel.background=element_rect(fill="NA", colour="black", size = 2, linetype='solid'),
          panel.grid.major=element_line(size=0.05,linetype='solid', colour="gray90"), 
          panel.grid.minor=element_line(size=0.05,linetype='solid', colour="gray90"),
          axis.text.y = element_text(size = 6),
          strip.text.x = element_text(size = 7),
          axis.title.y=element_blank()) +
    expand_limits(x=c(0, 5)) +
    scale_x_continuous(breaks=c(0,1,2,3,4,5)))
#legend.position = "none"))

png("Enrichment_genes_in_EAB_treatemnt only.png", res = 300, width = 2500, height = 1600)
plot (q6)# Make plot
dev.off() 

svg("Enrichment_B9VB20.svg", width = 7, height = 6)
plot (q6)# Make plot
dev.off() 
