library(ggplot2)
library(Matrix)
library(Seurat)
library(dplyr)
library(reshape2)
library(pheatmap)

norm_cross_col <- function(vec_value){
	vec_value_return <- (vec_value-mean(vec_value))/sd(vec_value)
}

MarkerInfo_HMDM <- read.delim("/public/ZhangJinLab/project_erv/DataForZhangLi/Cleandata/ZJU-ZL/singlecellfile/HMDM_M1M2_markers.txt", header = TRUE, stringsAsFactors = FALSE)
MarkerInfo_IPSDM <- read.delim("/public/ZhangJinLab/project_erv/DataForZhangLi/Cleandata/ZJU-ZL/singlecellfile/IPSDM_M1M2_markers.txt", header = TRUE, stringsAsFactors = FALSE)

MarkerInfo <- cbind(MarkerInfo_HMDM,MarkerInfo_IPSDM[,seq(3,12)])
M1M2Expr <- t(apply(MarkerInfo[,seq(3,14)],1,norm_cross_col))
M1M2order <- order(M1M2Expr[,1],decreasing=TRUE)
geneNames <- MarkerInfo[M1M2order,2]
M1M2Expr <- M1M2Expr[M1M2order,]
M1M2_plotdata_value <- as.vector(as.numeric(M1M2Expr))
M1M2_plotdata_genelabels <- rep(geneNames,12)
M1M2_plotdata_clulabels <- rep(c("HM_M1","HM_M2","IPS_M1_rep1","IPS_M1_rep2","IPS_M1_rep3","IPS_M1_rep4","IPS_M1_rep5","IPS_M2_rep1","IPS_M2_rep2","IPS_M2_rep3","IPS_M2_rep4","IPS_M2_rep5"),each=length(geneNames))

load("/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/AA/10X_LAH_InHouseAA.RData")
dataMatAA <- as.matrix(GetAssayData(object = mydata_total[["RNA"]], slot = "data"))
dataMatAA <- rowSums(dataMatAA)


load("/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/BB/10X_LAH_InHouseBB.RData")
dataMatBB <- as.matrix(GetAssayData(object = mydata_total[["RNA"]], slot = "data"))
dataMatBB <- rowSums(dataMatBB)


load("/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/DD/10X_LAH_InHouseDD.RData")
dataMatDD <- as.matrix(GetAssayData(object = mydata_total[["RNA"]], slot = "data"))
dataMatDD <- rowSums(dataMatDD)


AABB <- intersect(names(dataMatAA),names(dataMatBB))
AABBDD <- intersect(AABB,names(dataMatDD))

dataMat <- cbind(dataMatAA[AABBDD],dataMatBB[AABBDD],dataMatDD[AABBDD])
sampInfo <- factor(c("AA","BB","DD"),levels=c("AA","BB","DD"))

matchIndexes <- match(geneNames,rownames(dataMat))
dataMat <- as.data.frame(t(dataMat[matchIndexes,]))
plotdata_sc_value <- as.numeric(unlist(lapply(by(dataMat,sampInfo,colMeans),norm_cross_col)))
plotdata_sc_value <- c(M1M2_plotdata_value,plotdata_sc_value)
plotdata_sc_clulabels <- rep(c("AA","BB","DD"),each=ncol(dataMat))
plotdata_sc_clulabels <- c(M1M2_plotdata_clulabels,plotdata_sc_clulabels)
plotdata_sc_genelabels <- rep(geneNames,3)
plotdata_sc_genelabels <- c(M1M2_plotdata_genelabels,plotdata_sc_genelabels)
plotdata_sc <- as.data.frame(cbind(plotdata_sc_clulabels,plotdata_sc_genelabels,plotdata_sc_value))
plotdata_sc$plotdata_sc_genelabels <- factor(plotdata_sc_genelabels,levels=geneNames)
plotdata_sc$plotdata_sc_clulabels <- factor(plotdata_sc_clulabels,levels=c("HM_M1","IPS_M1_rep1","IPS_M1_rep2","IPS_M1_rep3","IPS_M1_rep4","IPS_M1_rep5","AA","BB","DD","HM_M2","IPS_M2_rep1","IPS_M2_rep2","IPS_M2_rep3","IPS_M2_rep4","IPS_M2_rep5"))
plotdata_sc$plotdata_sc_value <- as.numeric(plotdata_sc_value)

plotdata_sc <- acast(plotdata_sc, plotdata_sc_genelabels~plotdata_sc_clulabels, value.var="plotdata_sc_value")
plotdata_sc <- plotdata_sc[c(seq(1,29),32,35,seq(30,31),seq(33,34),seq(36,43)),]
genenames <- c("IL1A","IL1B","IL1RN","IL32","CD48","CD82","IL32A","CD80")
matchIndexes <- match(rownames(plotdata_sc),genenames)
plotdata_sc <- plotdata_sc[which(is.na(matchIndexes)),]

pdf(file="/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/M1M2_traditional_cluster_total.pdf",width=10,height=10)
pheatmap(plotdata_sc,scale="column",cluster_rows=FALSE,cluster_cols=TRUE,color = colorRampPalette(c("green","white","red"))(50),fontsize=4)
dev.off()

pdf(file="/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/M1M2_traditional_nocluster_total.pdf",width=10,height=10)
pheatmap(plotdata_sc,scale="column",cluster_rows=FALSE,cluster_cols=FALSE,color = colorRampPalette(c("green","white","red"))(50),fontsize=4)
dev.off()