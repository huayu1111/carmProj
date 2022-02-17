library(tidyverse)
library(ggplot2)

tpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

rawDataAA <- read.table(file="/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/AA/10X_LAH_InHouse_rawcounts.txt",row.names=1,header=TRUE)
CD86AA <- length(which(rawDataAA["CD86",] > 0))/ncol(rawDataAA)
CD68AA <- length(which(rawDataAA["CD68",] > 0))/ncol(rawDataAA)
CD80AA <- length(which(rawDataAA["CD80",] > 0))/ncol(rawDataAA)
CD206AA <- length(which(rawDataAA["MRC1",] > 0))/ncol(rawDataAA)
rawDataBB <- read.table(file="/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/BB/10X_LAH_InHouse_rawcounts.txt",row.names=1,header=TRUE)
CD86BB <- length(which(rawDataBB["CD86",] > 0))/ncol(rawDataBB)
CD68BB <- length(which(rawDataBB["CD68",] > 0))/ncol(rawDataBB)
CD80BB <- length(which(rawDataBB["CD80",] > 0))/ncol(rawDataBB)
CD206BB <- length(which(rawDataBB["MRC1",] > 0))/ncol(rawDataBB)
rawDataDD <- read.table(file="/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/DD/10X_LAH_InHouse_rawcounts.txt",row.names=1,header=TRUE)
CD86DD <- length(which(rawDataDD["CD86",] > 0))/ncol(rawDataDD)
CD68DD <- length(which(rawDataDD["CD68",] > 0))/ncol(rawDataDD)
CD80DD <- length(which(rawDataDD["CD80",] > 0))/ncol(rawDataDD)
CD206DD <- length(which(rawDataDD["MRC1",] > 0))/ncol(rawDataDD)

cond <- rep(c("AA","BB","DD"),each=4)
gene <- rep(c("CD86","CD68","CD80","CD206"),3)
perc <-c(CD86AA,CD68AA,CD80AA,CD206AA,CD86BB,CD68BB,CD80BB,CD206BB,CD86DD,CD68DD,CD80DD,CD206DD)

plotdata <- as.data.frame(cbind(cond,gene,perc),stringsAsFactors=FALSE)
colnames(plotdata) <- c("cond","gene","perc")
plotdata$cond <- factor(plotdata$cond,levels=c("AA","BB","DD"))
plotdata$gene <- factor(plotdata$gene,levels=c("CD86","CD80","CD68","CD206"))
plotdata$perc <- as.numeric(plotdata$perc)

pdf("/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/AA/M1M2marker_perc_add_CD68.pdf",width=7,height=5) 
p <- ggplot(data=plotdata, aes(x=gene, y=perc, fill=cond)) + geom_bar(stat="identity", color="black", position=position_dodge()) + scale_fill_manual(values=c("darkorchid3","goldenrod1","dodgerblue4"))
p <- p + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x=element_text(size=12), axis.title.y=element_text(size=12),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11))
print(p)
dev.off()


rawDataAA <- read.table(file="/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/AA/10X_LAH_InHouse_rawcounts.txt",row.names=1,header=TRUE)
rawDataBB <- read.table(file="/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/BB/10X_LAH_InHouse_rawcounts.txt",row.names=1,header=TRUE)
rawDataDD <- read.table(file="/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/DD/10X_LAH_InHouse_rawcounts.txt",row.names=1,header=TRUE)
expDataAA <- rowSums(rawDataAA)
expDataBB <- rowSums(rawDataBB)
expDataDD <- rowSums(rawDataDD)
matchIndexes <- match(rownames(rawDataAA),rownames(rawDataBB))
expDataAABB <- cbind(expDataAA[which(!is.na(matchIndexes))],expDataBB[matchIndexes[which(!is.na(matchIndexes))]])
matchIndexes <- match(rownames(expDataAABB),rownames(rawDataDD))
expData <- cbind(expDataAABB[which(!is.na(matchIndexes)),],expDataDD[matchIndexes[which(!is.na(matchIndexes))]])

colnames(expData) <- c("AA","BB","DD")
genemapData <- read.table(file="/public/ZhangJinLab/project_erv/refAnno/gencode.v28.genetype.tab",sep="\t",stringsAsFactors=FALSE,header=FALSE)
genelenData <- read.table(file="/public/ZhangJinLab/project_erv/DataForZhangLi/gene_cds_length.txt",sep=" ",stringsAsFactors=FALSE,header=FALSE)
matchIndexes <- match(rownames(expData),genelenData[,1])
expData <- expData[which(!is.na(matchIndexes)),]

expData <- tpm(expData,as.numeric(genelenData[matchIndexes[which(!is.na(matchIndexes))],2]))
matchIndexes <- match(rownames(expData),genemapData[,1])
expData <- cbind(genemapData[matchIndexes,3],expData)
write.table(expData,file="/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/AABBDD_gene_TPM.txt",sep='\t',quote=FALSE,row.names=T,col.names=T)


rawDataAA <- read.table(file="/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/AA/10X_LAH_InHouse_rawcounts.txt",row.names=1,header=TRUE)
rawDataBB <- read.table(file="/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/BB/10X_LAH_InHouse_rawcounts.txt",row.names=1,header=TRUE)
rawDataDD <- read.table(file="/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/DD/10X_LAH_InHouse_rawcounts.txt",row.names=1,header=TRUE)
expDataAACD86 <- rowSums(rawDataAA[,which(rawDataAA["CD86",] > 0)])
expDataAANonCD86 <- rowSums(rawDataAA[,which(rawDataAA["CD86",] == 0)])
expDataBBCD86 <- rowSums(rawDataBB[,which(rawDataBB["CD86",] > 0)])
expDataBBNonCD86 <- rowSums(rawDataBB[,which(rawDataBB["CD86",] == 0)])
expDataDDCD86 <- rowSums(rawDataDD[,which(rawDataDD["CD86",] > 0)])
expDataDDNonCD86 <- rowSums(rawDataDD[,which(rawDataDD["CD86",] == 0)])
matchIndexes <- match(rownames(rawDataAA),rownames(rawDataBB))
expDataAABB <- cbind(expDataAACD86[which(!is.na(matchIndexes))],expDataAANonCD86[which(!is.na(matchIndexes))],expDataBBCD86[matchIndexes[which(!is.na(matchIndexes))]],expDataBBNonCD86[matchIndexes[which(!is.na(matchIndexes))]])
matchIndexes <- match(rownames(expDataAABB),rownames(rawDataDD))
expData <- cbind(expDataAABB[which(!is.na(matchIndexes)),],expDataDDCD86[matchIndexes[which(!is.na(matchIndexes))]],expDataDDNonCD86[matchIndexes[which(!is.na(matchIndexes))]])

colnames(expData) <- c("AACD86","AANonCD86","BBCD86","BBNonCD86","DDCD86","DDNonCD86")
genemapData <- read.table(file="/public/ZhangJinLab/project_erv/refAnno/gencode.v28.genetype.tab",sep="\t",stringsAsFactors=FALSE,header=FALSE)
genelenData <- read.table(file="/public/ZhangJinLab/project_erv/DataForZhangLi/gene_cds_length.txt",sep=" ",stringsAsFactors=FALSE,header=FALSE)
matchIndexes <- match(rownames(expData),genemapData[,3])
expData <- expData[which(!is.na(matchIndexes)),]
rownames(expData) <- genemapData[matchIndexes[which(!is.na(matchIndexes))],1]
matchIndexes <- match(rownames(expData),genelenData[,1])
expData <- expData[which(!is.na(matchIndexes)),]

expData <- tpm(expData,as.numeric(genelenData[matchIndexes[which(!is.na(matchIndexes))],2]))
matchIndexes <- match(rownames(expData),genemapData[,1])
expData <- cbind(genemapData[matchIndexes,3],expData)
write.table(expData,file="/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/AABBDDCD86_gene_TPM.txt",sep='\t',quote=FALSE,row.names=T,col.names=T)


expDataAACD80 <- rowSums(rawDataAA[,which(rawDataAA["CD80",] > 0)])
expDataAANonCD80 <- rowSums(rawDataAA[,which(rawDataAA["CD80",] == 0)])
expDataBBCD80 <- rowSums(rawDataBB[,which(rawDataBB["CD80",] > 0)])
expDataBBNonCD80 <- rowSums(rawDataBB[,which(rawDataBB["CD80",] == 0)])
expDataDDCD80 <- rowSums(rawDataDD[,which(rawDataDD["CD80",] > 0)])
expDataDDNonCD80 <- rowSums(rawDataDD[,which(rawDataDD["CD80",] == 0)])
matchIndexes <- match(rownames(rawDataAA),rownames(rawDataBB))
expDataAABB <- cbind(expDataAACD80[which(!is.na(matchIndexes))],expDataAANonCD80[which(!is.na(matchIndexes))],expDataBBCD80[matchIndexes[which(!is.na(matchIndexes))]],expDataBBNonCD80[matchIndexes[which(!is.na(matchIndexes))]])
matchIndexes <- match(rownames(expDataAABB),rownames(rawDataDD))
expData <- cbind(expDataAABB[which(!is.na(matchIndexes)),],expDataDDCD80[matchIndexes[which(!is.na(matchIndexes))]],expDataDDNonCD80[matchIndexes[which(!is.na(matchIndexes))]])

colnames(expData) <- c("AACD80","AANonCD80","BBCD80","BBNonCD80","DDCD80","DDNonCD80")
genemapData <- read.table(file="/public/ZhangJinLab/project_erv/refAnno/gencode.v28.genetype.tab",sep="\t",stringsAsFactors=FALSE,header=FALSE)
genelenData <- read.table(file="/public/ZhangJinLab/project_erv/DataForZhangLi/gene_cds_length.txt",sep=" ",stringsAsFactors=FALSE,header=FALSE)
matchIndexes <- match(rownames(expData),genemapData[,3])
expData <- expData[which(!is.na(matchIndexes)),]
rownames(expData) <- genemapData[matchIndexes[which(!is.na(matchIndexes))],1]
matchIndexes <- match(rownames(expData),genelenData[,1])
expData <- expData[which(!is.na(matchIndexes)),]

expData <- tpm(expData,as.numeric(genelenData[matchIndexes[which(!is.na(matchIndexes))],2]))
matchIndexes <- match(rownames(expData),genemapData[,1])
expData <- cbind(genemapData[matchIndexes,3],expData)
write.table(expData,file="/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/AABBDDCD80_gene_TPM.txt",sep='\t',quote=FALSE,row.names=T,col.names=T)


expDataAACD68 <- rowSums(rawDataAA[,which(rawDataAA["CD68",] > 0)])
expDataAANonCD68 <- rowSums(rawDataAA[,which(rawDataAA["CD68",] == 0)])
expDataBBCD68 <- rowSums(rawDataBB[,which(rawDataBB["CD68",] > 0)])
expDataBBNonCD68 <- rowSums(rawDataBB[,which(rawDataBB["CD68",] == 0)])
expDataDDCD68 <- rowSums(rawDataDD[,which(rawDataDD["CD68",] > 0)])
expDataDDNonCD68 <- rowSums(rawDataDD[,which(rawDataDD["CD68",] == 0)])
matchIndexes <- match(rownames(rawDataAA),rownames(rawDataBB))
expDataAABB <- cbind(expDataAACD68[which(!is.na(matchIndexes))],expDataAANonCD68[which(!is.na(matchIndexes))],expDataBBCD68[matchIndexes[which(!is.na(matchIndexes))]],expDataBBNonCD68[matchIndexes[which(!is.na(matchIndexes))]])
matchIndexes <- match(rownames(expDataAABB),rownames(rawDataDD))
expData <- cbind(expDataAABB[which(!is.na(matchIndexes)),],expDataDDCD68[matchIndexes[which(!is.na(matchIndexes))]],expDataDDNonCD68[matchIndexes[which(!is.na(matchIndexes))]])

colnames(expData) <- c("AACD68","AANonCD68","BBCD68","BBNonCD68","DDCD68","DDNonCD68")
genemapData <- read.table(file="/public/ZhangJinLab/project_erv/refAnno/gencode.v28.genetype.tab",sep="\t",stringsAsFactors=FALSE,header=FALSE)
genelenData <- read.table(file="/public/ZhangJinLab/project_erv/DataForZhangLi/gene_cds_length.txt",sep=" ",stringsAsFactors=FALSE,header=FALSE)
matchIndexes <- match(rownames(expData),genemapData[,3])
expData <- expData[which(!is.na(matchIndexes)),]
rownames(expData) <- genemapData[matchIndexes[which(!is.na(matchIndexes))],1]
matchIndexes <- match(rownames(expData),genelenData[,1])
expData <- expData[which(!is.na(matchIndexes)),]

expData <- tpm(expData,as.numeric(genelenData[matchIndexes[which(!is.na(matchIndexes))],2]))
matchIndexes <- match(rownames(expData),genemapData[,1])
expData <- cbind(genemapData[matchIndexes,3],expData)
write.table(expData,file="/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/AABBDDCD68_gene_TPM.txt",sep='\t',quote=FALSE,row.names=T,col.names=T)


expDataAAMRC1 <- rowSums(rawDataAA[,which(rawDataAA["MRC1",] > 0)])
expDataAANonMRC1 <- rowSums(rawDataAA[,which(rawDataAA["MRC1",] == 0)])
expDataBBMRC1 <- rowSums(rawDataBB[,which(rawDataBB["MRC1",] > 0)])
expDataBBNonMRC1 <- rowSums(rawDataBB[,which(rawDataBB["MRC1",] == 0)])
expDataDDMRC1 <- rowSums(rawDataDD[,which(rawDataDD["MRC1",] > 0)])
expDataDDNonMRC1 <- rowSums(rawDataDD[,which(rawDataDD["MRC1",] == 0)])
matchIndexes <- match(rownames(rawDataAA),rownames(rawDataBB))
expDataAABB <- cbind(expDataAAMRC1[which(!is.na(matchIndexes))],expDataAANonMRC1[which(!is.na(matchIndexes))],expDataBBMRC1[matchIndexes[which(!is.na(matchIndexes))]],expDataBBNonMRC1[matchIndexes[which(!is.na(matchIndexes))]])
matchIndexes <- match(rownames(expDataAABB),rownames(rawDataDD))
expData <- cbind(expDataAABB[which(!is.na(matchIndexes)),],expDataDDMRC1[matchIndexes[which(!is.na(matchIndexes))]],expDataDDNonMRC1[matchIndexes[which(!is.na(matchIndexes))]])

colnames(expData) <- c("AAMRC1","AANonMRC1","BBMRC1","BBNonMRC1","DDMRC1","DDNonMRC1")
genemapData <- read.table(file="/public/ZhangJinLab/project_erv/refAnno/gencode.v28.genetype.tab",sep="\t",stringsAsFactors=FALSE,header=FALSE)
genelenData <- read.table(file="/public/ZhangJinLab/project_erv/DataForZhangLi/gene_cds_length.txt",sep=" ",stringsAsFactors=FALSE,header=FALSE)
matchIndexes <- match(rownames(expData),genemapData[,3])
expData <- expData[which(!is.na(matchIndexes)),]
rownames(expData) <- genemapData[matchIndexes[which(!is.na(matchIndexes))],1]
matchIndexes <- match(rownames(expData),genelenData[,1])
expData <- expData[which(!is.na(matchIndexes)),]

expData <- tpm(expData,as.numeric(genelenData[matchIndexes[which(!is.na(matchIndexes))],2]))
matchIndexes <- match(rownames(expData),genemapData[,1])
expData <- cbind(genemapData[matchIndexes,3],expData)
write.table(expData,file="/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/AABBDDCD206_gene_TPM.txt",sep='\t',quote=FALSE,row.names=T,col.names=T)       