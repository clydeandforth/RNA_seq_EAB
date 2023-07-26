library(edgeR)

setwd("/Users/krg114/Downloads/RNA_seq_EAB/EAB_ADB_Shared_HISAT2_featureCount/")


## Tree B3

EdgeRFile1 <- "b3_counts.csv"
countData_b3 <- read.table(EdgeRFile1, header = T, sep = ";", row.names = 1)
head(countData_b3)
str(countData_b3)

mygroups <- c("C1","C1","C1","C2","C2","C2")
y <- DGEList(counts=countData_b3, genes=rownames(countData_b3), group = mygroups)
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
result_edgeR <- as.data.frame(topTags(et, n=nrow(countData_b3)))

table(result_edgeR$FDR < 0.05)

plot(result_edgeR$logFC, -log10(result_edgeR$FDR), col=ifelse(result_edgeR$FDR<0.05,"red","black"),main="FDR volcano plot",xlab="log2FC",ylab="-log10(FDR)")
hist(result_edgeR$PValue, breaks=20, xlab="P-Value", col="royalblue", ylab="Frequency", main="P-value distribution")


comp_table <- merge(as.data.frame(res1), result_edgeR, by="row.names")

table("DESeq2" = res1$padj < 0.05, "edgeR" = comp_table$FDR < 0.05)
data_b3<-as.data.frame(table("DESeq2" = res1$padj < 0.05, "edgeR" = comp_table$FDR < 0.05))
total<-sum(data_b3$Freq)   
add_b3<-data_b3$Freq[1]+data_b3$Freq[4]
100/total*add_b3
# B3 = 61 % concordance
#########################
EdgeRFile1 <- "b8_counts.csv"
countData_b8 <- read.table(EdgeRFile1, header = T, sep = ",", row.names = 1)
head(countData_b8)
str(countData_b8)

mygroups <- c("C1","C1","C1","C2","C2","C2")
y <- DGEList(counts=countData_b8, genes=rownames(countData_b8), group = mygroups)
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
result_edgeR <- as.data.frame(topTags(et, n=nrow(countData_b8)))

table(result_edgeR$FDR < 0.05)

plot(result_edgeR$logFC, -log10(result_edgeR$FDR), col=ifelse(result_edgeR$FDR<0.05,"red","black"),main="FDR volcano plot",xlab="log2FC",ylab="-log10(FDR)")
hist(result_edgeR$PValue, breaks=20, xlab="P-Value", col="royalblue", ylab="Frequency", main="P-value distribution")

#Compare B3 result in DeSeq2 and EdgeR
comp_table <- merge(as.data.frame(res1), result_edgeR, by="row.names")

table("DESeq2" = res1$padj < 0.05, "edgeR" = comp_table$FDR < 0.05)
data_b8<-as.data.frame(table("DESeq2" = res1$padj < 0.05, "edgeR" = comp_table$FDR < 0.05))
total<-sum(data_b8$Freq)    

add_b8<-data_b8$Freq[1]+data_b8$Freq[4]
100/total*add_b8
# B8 = 60 % concordance
##############################

#########################
EdgeRFile1 <- "b9_counts.csv"
countData_b9 <- read.table(EdgeRFile1, header = T, sep = ",", row.names = 1)
head(countData_b9)
str(countData_b9)

mygroups <- c("C1","C1","C1","C2","C2","C2")
y <- DGEList(counts=countData_b9, genes=rownames(countData_b9), group = mygroups)
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
result_edgeR <- as.data.frame(topTags(et, n=nrow(countData_b9)))

table(result_edgeR$FDR < 0.05)

plot(result_edgeR$logFC, -log10(result_edgeR$FDR), col=ifelse(result_edgeR$FDR<0.05,"red","black"),main="FDR volcano plot",xlab="log2FC",ylab="-log10(FDR)")
hist(result_edgeR$PValue, breaks=20, xlab="P-Value", col="royalblue", ylab="Frequency", main="P-value distribution")

#Compare B3 result in DeSeq2 and EdgeR
comp_table <- merge(as.data.frame(res1), result_edgeR, by="row.names")

table("DESeq2" = res1$padj < 0.05, "edgeR" = comp_table$FDR < 0.05)
data_b9<-as.data.frame(table("DESeq2" = res1$padj < 0.05, "edgeR" = comp_table$FDR < 0.05))
total<-sum(data_b9$Freq)    

add_b9<-data_b9$Freq[1]+data_b9$Freq[4]
100/total*add_b9
# B9 = 60 % concordance

#########################
EdgeRFile1 <- "b20_counts.csv"
countData_b20 <- read.table(EdgeRFile1, header = T, sep = ",", row.names = 1)
head(countData_b20)
str(countData_b20)

mygroups <- c("C1","C1","C1","C2","C2","C2")
y <- DGEList(counts=countData_b20, genes=rownames(countData_b20), group = mygroups)
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
result_edgeR <- as.data.frame(topTags(et, n=nrow(countData_b20)))

table(result_edgeR$FDR < 0.05)

plot(result_edgeR$logFC, -log10(result_edgeR$FDR), col=ifelse(result_edgeR$FDR<0.05,"red","black"),main="FDR volcano plot",xlab="log2FC",ylab="-log10(FDR)")
hist(result_edgeR$PValue, breaks=20, xlab="P-Value", col="royalblue", ylab="Frequency", main="P-value distribution")

#Compare B3 result in DeSeq2 and EdgeR
comp_table <- merge(as.data.frame(res1), result_edgeR, by="row.names")

table("DESeq2" = res1$padj < 0.05, "edgeR" = comp_table$FDR < 0.05)
data_b20<-as.data.frame(table("DESeq2" = res1$padj < 0.05, "edgeR" = comp_table$FDR < 0.05))
total<-sum(data_b20$Freq)    

add_b20<-data_b20$Freq[1]+data_b20$Freq[4]
100/total*add_b20
# B20 = 80 % concordance

