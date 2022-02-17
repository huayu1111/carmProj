library(Matrix)
library(Seurat)
library(dplyr)
library(reshape2)

norm_cross_col <- function(vec_value){
	vec_value_return <- (vec_value-mean(vec_value))/sd(vec_value)
}


load("/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/AA/10X_LAH_InHouseAA.RData")
dataMatAA <- as.matrix(GetAssayData(object = mydata_total[["RNA"]], slot = "data"))
dataMatAACD86 <- rowSums(dataMatAA[,which(dataMatAA["CD86",] > 0)])
dataMatAANonCD86 <- rowSums(dataMatAA[,which(dataMatAA["CD86",] == 0)])

load("/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/BB/10X_LAH_InHouseBB.RData")
dataMatBB <- as.matrix(GetAssayData(object = mydata_total[["RNA"]], slot = "data"))
dataMatBBCD86 <- rowSums(dataMatBB[,which(dataMatBB["CD86",] > 0)])
dataMatBBNonCD86 <- rowSums(dataMatBB[,which(dataMatBB["CD86",] == 0)])

load("/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/DD/10X_LAH_InHouseDD.RData")
dataMatDD <- as.matrix(GetAssayData(object = mydata_total[["RNA"]], slot = "data"))
dataMatDDCD86 <- rowSums(dataMatDD[,which(dataMatDD["CD86",] > 0)])
dataMatDDNonCD86 <- rowSums(dataMatDD[,which(dataMatDD["CD86",] == 0)])

AABB <- intersect(names(dataMatAACD86),names(dataMatBBCD86))
AABBDD <- intersect(AABB,names(dataMatDDCD86))

sc_expdata <- cbind(dataMatAACD86[AABBDD],dataMatAANonCD86[AABBDD],dataMatBBCD86[AABBDD],dataMatBBNonCD86[AABBDD],dataMatDDCD86[AABBDD],dataMatDDNonCD86[AABBDD])

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
colnames(plotData) <- c("AACD86","AANonCD86","BBCD86","BBNonCD86","DDCD86","DDNonCD86")
plotData <- apply(plotData,2,norm_cross_col)


library(ggplot2)
plotData <- melt(plotData)
colnames(plotData) <- c("Pathways","Clusters","Expression")
plotData$Pathways <- factor(plotData$Pathways,levels=geneClasses)
plotData$Clusters <- factor(plotData$Clusters,levels=c("AACD86","AANonCD86","BBCD86","BBNonCD86","DDCD86","DDNonCD86"))
p <- ggplot(data = plotData, aes(x = Clusters, y = Pathways)) + geom_tile(aes(fill = Expression), colour = "white") + scale_fill_gradient2(low="blue", high="red", guide="colorbar")
p <- p + theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),axis.text.x=element_text(size=16),axis.text.y=element_text(size=16),legend.title=element_text(size=14),legend.text=element_text(size=11))
p <- p + theme(axis.text.x = element_text(vjust = 0.6, angle = 45))
ggsave(file=paste("/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/heatmap_Clusters_Pathways_sub_CD86.pdf", sep=""),width=25,height=15)
