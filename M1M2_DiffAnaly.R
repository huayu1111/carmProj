library(DESeq2)
mapInfo <- read.table(file = "/public/ZhangJinLab/project_erv/DataForZhangLi/SRP039361/stringtiefile/gene_TPM.txt", sep = "\t", header = T, row.names = 1)
database_all <- read.table(file = "/public/ZhangJinLab/project_erv/DataForZhangLi/SRP039361/stringtiefile/gene_count_matrix.csv", sep = ",", header = T)
expdata <- database_all[,2:ncol(database_all)]
rownames(expdata) <- database_all[,1]
conditions <- factor(c(rep("HMDM_M1",1),rep("HMDM_M2",1),rep("HMDM_M1",5),rep("HMDM_M2",5)))
coldata <- data.frame(row.names = colnames(expdata[,seq(1,12)]),conditions)
dds <- DESeqDataSetFromMatrix(countData=expdata[,seq(1,12)], colData=coldata, design=~conditions)
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj),]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
matchIndexes <- match(resdata[,1],mapInfo[,1])
rownames(resdata) <- resdata[,1]
resdata[,1] <- mapInfo[matchIndexes,2]
colnames(resdata)[1] <- c("gene_name")
resdata <- resdata[order(resdata$padj),]
# resdata <- resdata[seq(1,50),]
# write.csv(resdata, file = "/public/ZhangJinLab/project_erv/DataForZhangLi/SRP039361/stringtiefile/M1_vs_M2.csv",row.names=T,col.names=T)
resdata <- resdata[seq(51,55),]
write.csv(resdata, file = "/public/ZhangJinLab/project_erv/DataForZhangLi/SRP039361/stringtiefile/M1_vs_M2_51-55.csv",row.names=T,col.names=T)

# conditions <- factor(c(rep("IPSDM_M1",5),rep("IPSDM_M2",5)))
# coldata <- data.frame(row.names = colnames(expdata[,seq(3,12)]), conditions)
# dds <- DESeqDataSetFromMatrix(countData=expdata[,seq(3,12)], colData=coldata, design=~conditions)
# dds <- DESeq(dds)
# res <- results(dds)
# res <- res[order(res$padj),]
# resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
# matchIndexes <- match(resdata[,1],mapInfo[,1])
# rownames(resdata) <- resdata[,1]
# resdata[,1] <- mapInfo[matchIndexes,2]
# colnames(resdata)[1] <- c("gene_name")
# write.csv(resdata, file = "/public/ZhangJinLab/project_erv/DataForZhangLi/SRP039361/stringtiefile/IPSDM_M1_vs_IPSDM_M2.csv",row.names=T,col.names=T)