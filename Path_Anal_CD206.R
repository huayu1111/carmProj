library(Matrix)
library(Seurat)
library(dplyr)
library(reshape2)

norm_cross_col <- function(vec_value){
	vec_value_return <- (vec_value-mean(vec_value))/sd(vec_value)
}


load("/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/AA/10X_LAH_InHouseAA.RData")
dataMatAA <- as.matrix(GetAssayData(object = mydata_total[["RNA"]], slot = "data"))
dataMatAAMRC1 <- rowSums(dataMatAA[,which(dataMatAA["MRC1",] > 0)])
dataMatAANonMRC1 <- rowSums(dataMatAA[,which(dataMatAA["MRC1",] == 0)])

load("/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/BB/10X_LAH_InHouseBB.RData")
dataMatBB <- as.matrix(GetAssayData(object = mydata_total[["RNA"]], slot = "data"))
dataMatBBMRC1 <- rowSums(dataMatBB[,which(dataMatBB["MRC1",] > 0)])
dataMatBBNonMRC1 <- rowSums(dataMatBB[,which(dataMatBB["MRC1",] == 0)])

load("/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/DD/10X_LAH_InHouseDD.RData")
dataMatDD <- as.matrix(GetAssayData(object = mydata_total[["RNA"]], slot = "data"))
dataMatDDMRC1 <- rowSums(dataMatDD[,which(dataMatDD["MRC1",] > 0)])
dataMatDDNonMRC1 <- rowSums(dataMatDD[,which(dataMatDD["MRC1",] == 0)])

AABB <- intersect(names(dataMatAAMRC1),names(dataMatBBMRC1))
AABBDD <- intersect(AABB,names(dataMatDDMRC1))

sc_expdata <- cbind(dataMatAAMRC1[AABBDD],dataMatAANonMRC1[AABBDD],dataMatBBMRC1[AABBDD],dataMatBBNonMRC1[AABBDD],dataMatDDMRC1[AABBDD],dataMatDDNonMRC1[AABBDD])

geneInfo <- read.delim(file="/public/ZhangJinLab/project_metabolism/annofile/iMac_pathways_zl.txt",header=FALSE,stringsAsFactors=FALSE)
geneClasses <- unique(geneInfo[,4])

pathData <- c()
usedgeneClasses <- c()
for(geneclass in geneClasses){
    matchIndexes <- which(geneInfo[,4]==geneclass)
    geneIndexes <- match(geneInfo[matchIndexes,2],rownames(sc_expdata))
    geneIndexes <- geneIndexes[which(!is.na(geneIndexes))]
    if(length(geneIndexes) > 2){
        usedgeneClasses <- c(usedgeneClasses,geneclass)
        if(length(pathData) == 0){
            pathData <- apply(sc_expdata[geneIndexes,],2,mean)
        }else{
            pathData <- rbind(pathData,apply(sc_expdata[geneIndexes,],2,mean))
        }
    }
}

usedgeneClasses <- c()
for(geneclass in geneClasses){
    matchIndexes <- which(geneInfo[,4]==geneclass)
    geneIndexes <- match(geneInfo[matchIndexes,2],rownames(sc_expdata))
    geneIndexes <- geneIndexes[which(!is.na(geneIndexes))]
    if(length(geneIndexes) > 2){
        usedgeneClasses <- c(usedgeneClasses,geneclass)
    }
}
rownames(pathData) <- usedgeneClasses

plotData <- pathData
rownames(plotData) <- usedgeneClasses
colnames(plotData) <- c("AAMRC1","AANonMRC1","BBMRC1","BBNonMRC1","DDMRC1","DDNonMRC1")
plotData <- apply(plotData,2,norm_cross_col)


library(ggplot2)
plotData <- melt(plotData)
colnames(plotData) <- c("Pathways","Clusters","Expression")
plotData$Pathways <- factor(plotData$Pathways,levels=geneClasses)
plotData$Clusters <- factor(plotData$Clusters,levels=c("AAMRC1","AANonMRC1","BBMRC1","BBNonMRC1","DDMRC1","DDNonMRC1"))
p <- ggplot(data = plotData, aes(x = Clusters, y = Pathways)) + geom_tile(aes(fill = Expression), colour = "white") + scale_fill_gradient2(low="blue", high="red", guide="colorbar")
p <- p + theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),axis.text.x=element_text(size=16),axis.text.y=element_text(size=16),legend.title=element_text(size=14),legend.text=element_text(size=11))
p <- p + theme(axis.text.x = element_text(vjust = 0.6, angle = 45))
ggsave(file=paste("/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/heatmap_Clusters_Pathways_sub_CD206.pdf", sep=""),width=25,height=15)
