library(ggplot2)
library(sva)

maxvalue <- function(x){
    max_value <- max(x)
}

expValTrans <- function(expVec){
	expVec <- (expVec-mean(expVec))/sd(expVec)
	return(expVec)
}

norm_cross_col <- function(vec_value){
	vec_value_return <- (vec_value-mean(vec_value))/sd(vec_value)
}

MarkerInfo_HMDM <- read.delim("/public/ZhangJinLab/backdata/lahdata/cdata/HMDM_M1M2_markers.txt", header = TRUE, stringsAsFactors = FALSE)
MarkerInfo_IPSDM <- read.delim("/public/ZhangJinLab/backdata/lahdata/cdata/IPSDM_M1M2_markers.txt", header = TRUE, stringsAsFactors = FALSE)

MarkerInfo <- cbind(MarkerInfo_HMDM,MarkerInfo_IPSDM[,seq(3,12)])
M1M2Expr <- t(apply(MarkerInfo[,seq(3,14)],1,norm_cross_col))
M1M2order <- order(M1M2Expr[,1],decreasing=TRUE)
geneNames <- MarkerInfo[M1M2order,2]

expData_SRP103864 <- read.table(file="/public/ZhangJinLab/backdata/lahdata/cdata/SRP103864/gene_TPM.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE)
rownames(expData_SRP103864) <- expData_SRP103864[,1]
matchIndexes <- match(geneNames,expData_SRP103864[,2])
geneEnsembleIDs <- expData_SRP103864[matchIndexes,1]
geneEnsembleIDs <- gsub("\\.\\d+","",geneEnsembleIDs)

EnsembIDs1 <- read.table(file="/public/ZhangJinLab/backdata/lahdata/M1statgenes.txt",sep="\t",header=FALSE,stringsAsFactors=FALSE)
EnsembIDs2 <- read.table(file="/public/ZhangJinLab/backdata/lahdata/M2statgenes.txt",sep="\t",header=FALSE,stringsAsFactors=FALSE)
EnsembIDs <- rbind(EnsembIDs1,EnsembIDs2)
geneEnsembleIDs <- unique(c(geneEnsembleIDs,as.character(EnsembIDs[,2])))
geneNames <- expData_SRP103864[match(geneEnsembleIDs,gsub("\\.\\d+","",rownames(expData_SRP103864))),2]
rGNames <- c("IL6","CD300E","RETNLB")
matchIndexes <- setdiff(seq(1,length(geneNames)),match(rGNames,geneNames))
geneEnsembleIDs <- geneEnsembleIDs[matchIndexes]
geneNames <- geneNames[matchIndexes] 
 
sampInfo <- read.table(file="/public/ZhangJinLab/backdata/lahdata/cdata/SRP103864/sampInfo_SRP103864.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE)
sampInfo[,2] <- factor(sampInfo[,2],levels=c("Non-polarized macrophages rep.1","Non-polarized macrophages rep.2","Macrophages polarized with IFN-gamma and LPS rep.1","Macrophages polarized with IFN-gamma and LPS rep.2","Macrophages polarized with IL-4 and IL-13 rep.1","Macrophages polarized with IL-4 and IL-13 rep.2","Macrophages polarized with IL-10 rep.1","Macrophages polarized with IL-10 rep.2"))
sampInfo <- sampInfo[order(sampInfo[,2]),]
matchIndexes <- match(sampInfo[,1],colnames(expData_SRP103864))
expData_SRP103864 <- expData_SRP103864[,matchIndexes]
colnames(expData_SRP103864) <- as.character(sampInfo[,2])

expData_SRP039361 <- read.table(file="/public/ZhangJinLab/backdata/lahdata/cdata/SRP039361/gene_TPM.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE)
rownames(expData_SRP039361) <- expData_SRP039361[,1]
sampInfo <- read.table(file="/public/ZhangJinLab/backdata/lahdata/cdata/SRP039361/sampInfo_SRP039361.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE)
sampInfo[,2] <- factor(sampInfo[,2],levels=c("HMDM M1 rep1","HMDM M1 rep2","HMDM M2 rep1","HMDM M2 rep2","HMDM MAC rep1","HMDM MAC rep2","HMDM MAC rep3","IPSDM M1 rep1","IPSDM M1 rep2","IPSDM M1 rep3","IPSDM M1 rep4","IPSDM M1 rep5","IPSDM M1 rep6","IPSDM M2 rep1","IPSDM M2 rep3","IPSDM M2 rep4","IPSDM M2 rep5","IPSDM M2 rep6","IPSDM MAC rep1","IPSDM MAC rep2","IPSDM MAC rep3","IPSDM MAC rep4","iPS rep1","iPS rep2","iPS rep5","iPS rep6"))
sampInfo <- sampInfo[order(sampInfo[,2]),]
matchIndexes <- match(sampInfo[,1],colnames(expData_SRP039361))
matchIndexes_have <- matchIndexes[which(!is.na(matchIndexes))]
expData_SRP039361 <- expData_SRP039361[,matchIndexes_have]
colnames(expData_SRP039361) <- as.character(sampInfo[which(!is.na(matchIndexes)),2])

expData_InHouse <- read.table(file="/public/ZhangJinLab/backdata/lahdata/cdata/AABBDD_gene_TPM.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE)
allGeneNames <- expData_InHouse[,1]
expData_InHouse <- expData_InHouse[,seq(2,ncol(expData_InHouse))]

matchIndexesA <- match(rownames(expData_InHouse),gsub("\\.\\d+","",rownames(expData_SRP103864)))
matchIndexesB <- match(rownames(expData_InHouse),gsub("\\.\\d+","",rownames(expData_SRP039361)))

expData <- cbind(expData_InHouse,expData_SRP103864[matchIndexesA,],expData_SRP039361[matchIndexesB,])
expData <- expData[,c("AA","BB","DD","IPSDM M1 rep1","IPSDM M1 rep2","IPSDM M1 rep3","IPSDM M2 rep1","IPSDM M2 rep6","IPSDM MAC rep2","IPSDM MAC rep3","IPSDM MAC rep4")]
expData <- expData[which(apply(expData,1,maxvalue) > 1),]
batch <- as.factor(rep(c("InHouse","SRP039361"), c(3,8)))
modcombat = model.matrix(~1,data=batch)
combat_edata = ComBat(dat=expData, batch=batch, mod=modcombat, par.prior=TRUE)
# combat_edata <- apply(combat_edata,1,expValTrans)

matchIndexes <- match(geneEnsembleIDs,rownames(combat_edata))
combat_edata <- combat_edata[matchIndexes[which(!is.na(matchIndexes))],]
rownames(combat_edata) <- geneNames[which(!is.na(matchIndexes))]

stageNames <- rep(c("AA","BB","DD","IPSDM M1","IPSDM M2","IPSDM MAC"),c(1,1,1,3,2,3))
PCAFitClass <- prcomp(t(combat_edata))
max3PCs <- PCAFitClass$x[,1:3]
x <- (max3PCs[,1]-min(max3PCs[,1]))/(max(max3PCs[,1])-min(max3PCs[,1]))
y <- (max3PCs[,2]-min(max3PCs[,2]))/(max(max3PCs[,2])-min(max3PCs[,2]))
plotdata <- data.frame(x=x,y=y,pointColor=stageNames)

p <- ggplot(plotdata,aes(x=x,y=y,colour=factor(pointColor,levels=c("AA","BB","DD","IPSDM M1","IPSDM M2","IPSDM MAC")),shape=factor(pointColor,levels=c("AA","BB","DD","IPSDM M1","IPSDM M2","IPSDM MAC")))) + geom_point(size=3) + scale_shape_manual(values=c(1,2,3,4,5,6))
p <- p + xlab("PC1") + ylab("PC2")
p <- p + theme(legend.title=element_blank())
p <- p + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11))
ggsave("/public/ZhangJinLab/backdata/lahdata/cdata/PCA_2DRes_PC1_PC2_ALL.pdf",width=4,height=3)

x <- (max3PCs[,1]-min(max3PCs[,1]))/(max(max3PCs[,1])-min(max3PCs[,1]))
y <- (max3PCs[,3]-min(max3PCs[,3]))/(max(max3PCs[,3])-min(max3PCs[,3]))
plotdata <- data.frame(x=x,y=y,pointColor=stageNames)

p <- ggplot(plotdata,aes(x=x,y=y,colour=factor(pointColor,levels=c("AA","BB","DD","IPSDM M1","IPSDM M2","IPSDM MAC")),shape=factor(pointColor,levels=c("AA","BB","DD","IPSDM M1","IPSDM M2","IPSDM MAC")))) + geom_point(size=3) + scale_shape_manual(values=c(1,2,3,4,5,6))
p <- p + xlab("PC1") + ylab("PC3")
p <- p + theme(legend.title=element_blank())
p <- p + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11))
ggsave("/public/ZhangJinLab/backdata/lahdata/cdata/PCA_2DRes_PC1_PC3_ALL.pdf",width=4,height=3)

x <- (max3PCs[,2]-min(max3PCs[,2]))/(max(max3PCs[,2])-min(max3PCs[,2]))
y <- (max3PCs[,3]-min(max3PCs[,3]))/(max(max3PCs[,3])-min(max3PCs[,3]))
plotdata <- data.frame(x=x,y=y,pointColor=stageNames)

p <- ggplot(plotdata,aes(x=x,y=y,colour=factor(pointColor,levels=c("AA","BB","DD","IPSDM M1","IPSDM M2","IPSDM MAC")),shape=factor(pointColor,levels=c("AA","BB","DD","IPSDM M1","IPSDM M2","IPSDM MAC")))) + geom_point(size=3) + scale_shape_manual(values=c(1,2,3,4,5,6))
p <- p + xlab("PC2") + ylab("PC3")
p <- p + theme(legend.title=element_blank())
p <- p + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11))
ggsave("/public/ZhangJinLab/backdata/lahdata/cdata/PCA_2DRes_PC2_PC3_ALL.pdf",width=4,height=3)


library(scatterplot3d)
x <- (max3PCs[,1]-min(max3PCs[,1]))/(max(max3PCs[,1])-min(max3PCs[,1]))
y <- (max3PCs[,2]-min(max3PCs[,2]))/(max(max3PCs[,2])-min(max3PCs[,2]))
z <- (max3PCs[,3]-min(max3PCs[,3]))/(max(max3PCs[,3])-min(max3PCs[,3]))
pdf(file="/public/ZhangJinLab/backdata/lahdata/cdata/PCA_3DRes_ALL.pdf",width=18,height=9)
par(xpd=TRUE,mgp=c(2,0.8,0),oma=c(0.2,0.2,0.2,0.2),lwd=3)
scatterplot3d(x,y,z,color=rep(rainbow(6),c(1,1,1,3,2,3)),pch=rep(seq(1,6),c(1,1,1,3,2,3)),xlab="PC1",ylab="PC2",zlab="PC3",mar=c(4,4,4,28),cex.axis=2,cex.lab=3,cex.symbols=2.5)
legend("topright",c("AA","BB","DD","IPSDM M1","IPSDM M2","IPSDM MAC"),inset=c(-0.47,-0.06),pch=seq(1,6),col=rainbow(6),cex=2,y.intersp=1.0)
dev.off()
