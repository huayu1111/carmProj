library(pheatmap)
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
expData <- expData[which(apply(expData,1,maxvalue) > 1),]
batch <- as.factor(rep(c("InHouse","SRP103864","SRP039361"), c(ncol(expData_InHouse),ncol(expData_SRP103864),ncol(expData_SRP039361))))
modcombat = model.matrix(~1,data=batch)
combat_edata = ComBat(dat=expData, batch=batch, mod=modcombat, par.prior=TRUE)
combat_edata <- combat_edata[,c("AA","BB","DD","Non-polarized macrophages rep.2","Macrophages polarized with IFN-gamma and LPS rep.2","IPSDM M1 rep1","IPSDM M1 rep2","IPSDM M1 rep3","IPSDM M2 rep1","IPSDM M2 rep6","IPSDM MAC rep2","IPSDM MAC rep3","IPSDM MAC rep4","iPS rep1","iPS rep2","iPS rep5","iPS rep6")]
# combat_edata <- apply(combat_edata,1,expValTrans)

pdf(file="/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/Cluster_Res.pdf",width=10,height=15)
p <- pheatmap(combat_edata,show_rownames=FALSE,scale="row",cluster_rows=FALSE,cluster_cols=TRUE,fontsize=10)
plot(p$tree_col)
dev.off()
