library(ggplot2)
library(Matrix)
library(Seurat)
library(dplyr)
library(reshape2)
library(pheatmap)

norm_cross_col <- function(vec_value){
	vec_value_return <- (vec_value-mean(vec_value))/sd(vec_value)
}

mapInfo <- read.table(file = "/public/ZhangJinLab/project_erv/DataForZhangLi/SRP039361/stringtiefile/gene_TPM.txt", sep = "\t", header = T, row.names = 1)
MarkerInfo <- read.delim("/public/ZhangJinLab/project_erv/DataForZhangLi/Cleandata/ZJU-ZL/singlecellfile/diffanaly_markers.txt", header = TRUE, stringsAsFactors = FALSE)
matchIndexes <- match(MarkerInfo[,2],mapInfo[,2])
MarkerInfo[,seq(3,ncol(MarkerInfo))] <- mapInfo[matchIndexes,c("SRR1182376","SRR1182388","SRR2939146","SRR1182390","SRR2910670","SRR2910671","SRR1182378","SRR1182392","SRR1182394","SRR2939152","SRR2910672","SRR2910673")]
M1M2Expr <- t(apply(MarkerInfo[,seq(3,14)],1,norm_cross_col))
M1M2order <- order(M1M2Expr[,1],decreasing=TRUE)
geneNames <- MarkerInfo[M1M2order,2]
M1M2Expr <- M1M2Expr[M1M2order,]
M1M2_plotdata_value <- as.vector(as.numeric(M1M2Expr))
M1M2_plotdata_genelabels <- rep(geneNames,12)
M1M2_plotdata_clulabels <- rep(c("HM_M1","IPS_M1_rep1","IPS_M1_rep2","IPS_M1_rep3","IPS_M1_rep4","IPS_M1_rep5","HM_M2","IPS_M2_rep1","IPS_M2_rep2","IPS_M2_rep3","IPS_M2_rep4","IPS_M2_rep5"),each=length(geneNames))

load("/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/AA/10X_LAH_InHouseAA.RData")
dataMatAA <- as.matrix(GetAssayData(object = mydata_total[["RNA"]], slot = "data"))
dataMatAACD80 <- rowSums(dataMatAA[,which(dataMatAA["CD80",] > 0)])
dataMatAANonCD80 <- rowSums(dataMatAA[,which(dataMatAA["CD80",] == 0)])

load("/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/BB/10X_LAH_InHouseBB.RData")
dataMatBB <- as.matrix(GetAssayData(object = mydata_total[["RNA"]], slot = "data"))
dataMatBBCD80 <- rowSums(dataMatBB[,which(dataMatBB["CD80",] > 0)])
dataMatBBNonCD80 <- rowSums(dataMatBB[,which(dataMatBB["CD80",] == 0)])

load("/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/DD/10X_LAH_InHouseDD.RData")
dataMatDD <- as.matrix(GetAssayData(object = mydata_total[["RNA"]], slot = "data"))
dataMatDDCD80 <- rowSums(dataMatDD[,which(dataMatDD["CD80",] > 0)])
dataMatDDNonCD80 <- rowSums(dataMatDD[,which(dataMatDD["CD80",] == 0)])

AABB <- intersect(names(dataMatAACD80),names(dataMatBBCD80))
AABBDD <- intersect(AABB,names(dataMatDDCD80))

dataMat <- cbind(dataMatAACD80[AABBDD],dataMatAANonCD80[AABBDD],dataMatBBCD80[AABBDD],dataMatBBNonCD80[AABBDD],dataMatDDCD80[AABBDD],dataMatDDNonCD80[AABBDD])
sampInfo <- factor(c("AACD80","AANonCD80","BBCD80","BBNonCD80","DDCD80","DDNonCD80"),levels=c("AACD80","AANonCD80","BBCD80","BBNonCD80","DDCD80","DDNonCD80"))

matchIndexes <- match(geneNames,rownames(dataMat))
dataMat <- as.data.frame(t(dataMat[matchIndexes,]))
plotdata_sc_value <- as.numeric(unlist(lapply(by(dataMat,sampInfo,colMeans),norm_cross_col)))
plotdata_sc_value <- c(M1M2_plotdata_value,plotdata_sc_value)
plotdata_sc_clulabels <- rep(c("AACD80","AANonCD80","BBCD80","BBNonCD80","DDCD80","DDNonCD80"),each=ncol(dataMat))
plotdata_sc_clulabels <- c(M1M2_plotdata_clulabels,plotdata_sc_clulabels)
plotdata_sc_genelabels <- rep(geneNames,6)
plotdata_sc_genelabels <- c(M1M2_plotdata_genelabels,plotdata_sc_genelabels)
plotdata_sc <- as.data.frame(cbind(plotdata_sc_clulabels,plotdata_sc_genelabels,plotdata_sc_value))
plotdata_sc$plotdata_sc_genelabels <- factor(plotdata_sc_genelabels,levels=geneNames)
plotdata_sc$plotdata_sc_clulabels <- factor(plotdata_sc_clulabels,levels=c("HM_M1","IPS_M1_rep1","IPS_M1_rep2","IPS_M1_rep3","IPS_M1_rep4","IPS_M1_rep5","AACD80","AANonCD80","BBCD80","BBNonCD80","DDCD80","DDNonCD80","HM_M2","IPS_M2_rep1","IPS_M2_rep2","IPS_M2_rep3","IPS_M2_rep4","IPS_M2_rep5"))
plotdata_sc$plotdata_sc_value <- as.numeric(plotdata_sc_value)

plotdata_sc <- acast(plotdata_sc, plotdata_sc_genelabels~plotdata_sc_clulabels, value.var="plotdata_sc_value")

pdf(file="/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/M1M2_diffanaly_cluster_CD80.pdf",width=10,height=10)
pheatmap(plotdata_sc,scale="none",cluster_rows=FALSE,cluster_cols=TRUE,breaks = c(seq(-2.5,0,0.1),seq(0.1,2.5,0.1)),color = colorRampPalette(c("green","white","red"))(50),fontsize=4)
dev.off()

pdf(file="/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/M1M2_diffanaly_nocluster_CD80.pdf",width=10,height=10)
pheatmap(plotdata_sc,scale="none",cluster_rows=FALSE,cluster_cols=FALSE,breaks = c(seq(-2.5,0,0.1),seq(0.1,2.5,0.1)),color = colorRampPalette(c("green","white","red"))(50),fontsize=4)
dev.off()