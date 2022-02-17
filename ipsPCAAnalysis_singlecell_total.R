library(ggplot2)
library(scatterplot3d)
library(sva)

maxvalue <- function(x){
    max_value <- max(x)
}

expValTrans <- function(expVec){
	expVec <- (expVec-mean(expVec))/sd(expVec)
	return(expVec)
}

expData_SRP103864 <- read.table(file="/public/ZhangJinLab/project_erv/DataForZhangLi/SRP103864/stringtiefile/gene_TPM.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE)
rownames(expData_SRP103864) <- expData_SRP103864[,1]
sampInfo <- read.table(file="/public/ZhangJinLab/project_erv/DataForZhangLi/SRP103864/stringtiefile/sampInfo_SRP103864.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE)
sampInfo[,2] <- factor(sampInfo[,2],levels=c("Non-polarized macrophages rep.1","Non-polarized macrophages rep.2","Macrophages polarized with IFN-gamma and LPS rep.1","Macrophages polarized with IFN-gamma and LPS rep.2","Macrophages polarized with IL-4 and IL-13 rep.1","Macrophages polarized with IL-4 and IL-13 rep.2","Macrophages polarized with IL-10 rep.1","Macrophages polarized with IL-10 rep.2"))
sampInfo <- sampInfo[order(sampInfo[,2]),]
matchIndexes <- match(sampInfo[,1],colnames(expData_SRP103864))
expData_SRP103864 <- expData_SRP103864[,matchIndexes]
colnames(expData_SRP103864) <- as.character(sampInfo[,2])

expData_SRP039361 <- read.table(file="/public/ZhangJinLab/project_erv/DataForZhangLi/SRP039361/stringtiefile/gene_TPM.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE)
rownames(expData_SRP039361) <- expData_SRP039361[,1]
sampInfo <- read.table(file="/public/ZhangJinLab/project_erv/DataForZhangLi/SRP039361/stringtiefile/sampInfo_SRP039361.txt",header=FALSE,sep="\t",stringsAsFactors=FALSE)
sampInfo[,2] <- factor(sampInfo[,2],levels=c("HMDM M1 rep1","HMDM M1 rep2","HMDM M2 rep1","HMDM M2 rep2","HMDM MAC rep1","HMDM MAC rep2","HMDM MAC rep3","IPSDM M1 rep1","IPSDM M1 rep2","IPSDM M1 rep3","IPSDM M1 rep4","IPSDM M1 rep5","IPSDM M1 rep6","IPSDM M2 rep1","IPSDM M2 rep3","IPSDM M2 rep4","IPSDM M2 rep5","IPSDM M2 rep6","IPSDM MAC rep1","IPSDM MAC rep2","IPSDM MAC rep3","IPSDM MAC rep4","iPS rep1","iPS rep2","iPS rep5","iPS rep6"))
sampInfo <- sampInfo[order(sampInfo[,2]),]
matchIndexes <- match(sampInfo[,1],colnames(expData_SRP039361))
matchIndexes_have <- matchIndexes[which(!is.na(matchIndexes))]
expData_SRP039361 <- expData_SRP039361[,matchIndexes_have]
colnames(expData_SRP039361) <- as.character(sampInfo[which(!is.na(matchIndexes)),2])

expData_InHouse <- read.table(file="/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/AABBDD_gene_TPM.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE)
expData_InHouse <- expData_InHouse[,seq(2,ncol(expData_InHouse))]

matchIndexesA <- match(rownames(expData_InHouse),rownames(expData_SRP103864))
matchIndexesB <- match(rownames(expData_InHouse),rownames(expData_SRP039361))

expData <- cbind(expData_InHouse,expData_SRP103864[matchIndexesA,],expData_SRP039361[matchIndexesB,])
expData <- expData[which(apply(expData,1,maxvalue) > 0),]
batch <- as.factor(rep(c("InHouse","SRP103864","SRP039361"), c(ncol(expData_InHouse),ncol(expData_SRP103864),ncol(expData_SRP039361))))
modcombat = model.matrix(~1,data=batch)
combat_edata = ComBat(dat=expData, batch=batch, mod=modcombat, par.prior=TRUE)
combat_edata <- combat_edata[,c("AA","BB","DD","Non-polarized macrophages rep.2","Macrophages polarized with IFN-gamma and LPS rep.2","IPSDM M1 rep1","IPSDM M1 rep2","IPSDM M1 rep3","IPSDM M2 rep1","IPSDM M2 rep6","IPSDM MAC rep2","IPSDM MAC rep3","IPSDM MAC rep4","iPS rep1","iPS rep2","iPS rep5","iPS rep6")]
combat_edata <- apply(combat_edata,1,expValTrans)

stageNames <- rep(c("AA","BB","DD","Non-polarized macrophages","Macrophages polarized with IFN-gamma and LPS","IPSDM M1","IPSDM M2","IPSDM MAC","iPS"),c(1,1,1,1,1,3,2,3,4))
PCAFitClass <- prcomp(combat_edata)
max3PCs <- PCAFitClass$x[,1:3]
x <- (max3PCs[,1]-min(max3PCs[,1]))/(max(max3PCs[,1])-min(max3PCs[,1]))
y <- (max3PCs[,2]-min(max3PCs[,2]))/(max(max3PCs[,2])-min(max3PCs[,2]))
plotdata <- data.frame(x=x,y=y,pointColor=stageNames)

p <- ggplot(plotdata,aes(x=x,y=y,colour=factor(pointColor,levels=c("AA","BB","DD","Non-polarized macrophages","Macrophages polarized with IFN-gamma and LPS","IPSDM M1","IPSDM M2","IPSDM MAC","iPS")),shape=factor(pointColor,levels=c("AA","BB","DD","Non-polarized macrophages","Macrophages polarized with IFN-gamma and LPS","IPSDM M1","IPSDM M2","IPSDM MAC","iPS")))) + geom_point(size=3) + scale_shape_manual(values=c(1,2,3,4,5,6,7,8,9))
p <- p + xlab("PC1") + ylab("PC2")
p <- p + theme(legend.title=element_blank())
p <- p + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11))
ggsave("/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/PCA_2DRes_PC1_PC2.pdf",width=10,height=6)

x <- (max3PCs[,1]-min(max3PCs[,1]))/(max(max3PCs[,1])-min(max3PCs[,1]))
y <- (max3PCs[,3]-min(max3PCs[,3]))/(max(max3PCs[,3])-min(max3PCs[,3]))
plotdata <- data.frame(x=x,y=y,pointColor=stageNames)

p <- ggplot(plotdata,aes(x=x,y=y,colour=factor(pointColor,levels=c("AA","BB","DD","Non-polarized macrophages","Macrophages polarized with IFN-gamma and LPS","IPSDM M1","IPSDM M2","IPSDM MAC","iPS")),shape=factor(pointColor,levels=c("AA","BB","DD","Non-polarized macrophages","Macrophages polarized with IFN-gamma and LPS","IPSDM M1","IPSDM M2","IPSDM MAC","iPS")))) + geom_point(size=3) + scale_shape_manual(values=c(1,2,3,4,5,6,7,8,9))
p <- p + xlab("PC1") + ylab("PC3")
p <- p + theme(legend.title=element_blank())
p <- p + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11))
ggsave("/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/PCA_2DRes_PC1_PC3.pdf",width=10,height=6)

x <- (max3PCs[,2]-min(max3PCs[,2]))/(max(max3PCs[,2])-min(max3PCs[,2]))
y <- (max3PCs[,3]-min(max3PCs[,3]))/(max(max3PCs[,3])-min(max3PCs[,3]))
plotdata <- data.frame(x=x,y=y,pointColor=stageNames)

p <- ggplot(plotdata,aes(x=x,y=y,colour=factor(pointColor,levels=c("AA","BB","DD","Non-polarized macrophages","Macrophages polarized with IFN-gamma and LPS","IPSDM M1","IPSDM M2","IPSDM MAC","iPS")),shape=factor(pointColor,levels=c("AA","BB","DD","Non-polarized macrophages","Macrophages polarized with IFN-gamma and LPS","IPSDM M1","IPSDM M2","IPSDM MAC","iPS")))) + geom_point(size=3) + scale_shape_manual(values=c(1,2,3,4,5,6,7,8,9))
p <- p + xlab("PC2") + ylab("PC3")
p <- p + theme(legend.title=element_blank())
p <- p + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11))
ggsave("/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/PCA_2DRes_PC2_PC3.pdf",width=10,height=6)

x <- (max3PCs[,1]-min(max3PCs[,1]))/(max(max3PCs[,1])-min(max3PCs[,1]))
y <- (max3PCs[,2]-min(max3PCs[,2]))/(max(max3PCs[,2])-min(max3PCs[,2]))
z <- (max3PCs[,3]-min(max3PCs[,3]))/(max(max3PCs[,3])-min(max3PCs[,3]))
# 3D plot
# palette(rainbow(9))
# png(file="/public/ZhangJinLab/project_erv/DataForZhangLi/Cleandata/stringtiefile/PCA_3DRes.png",width=30,height=18,units="in",res=400)
pdf(file="/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/PCA_3DRes.pdf",width=18,height=9)
par(xpd=TRUE,mgp=c(2,0.8,0),oma=c(0.2,0.2,0.2,0.2),lwd=3)
scatterplot3d(x,y,z,color=rep(rainbow(9),c(1,1,1,1,1,3,2,3,4)),pch=rep(seq(1,9),c(1,1,1,1,1,3,2,3,4)),xlab="PC1",ylab="PC2",zlab="PC3",mar=c(4,4,4,28),cex.axis=2,cex.lab=3,cex.symbols=2.5)
legend("topright",c("AA","BB","DD","Non-polarized macrophages","Polarized with IFN-gamma and LPS","IPSDM M1","IPSDM M2","IPSDM MAC","iPS"),inset=c(-0.47,-0.01),pch=seq(1,9),col=rainbow(9),cex=2,y.intersp=1.0)
dev.off()