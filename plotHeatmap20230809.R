library(Seurat)
load("/public/ZhangJinLab/backdata/lahdata/AA/10X_LAH_InHouse.RData")
normDataAA <-  GetAssayData(object = mydata_total[["RNA"]], slot = "data")
load("/public/ZhangJinLab/backdata/lahdata/BB/10X_LAH_InHouse.RData")
normDataBB <-  GetAssayData(object = mydata_total[["RNA"]], slot = "data")
load("/public/ZhangJinLab/backdata/lahdata/DD/10X_LAH_InHouse.RData")
normDataDD <-  GetAssayData(object = mydata_total[["RNA"]], slot = "data")

# Part 1
EnsembIDs1 <- read.table(file="/public/ZhangJinLab/backdata/lahdata/M1statgenes.txt",sep="\t",header=FALSE,stringsAsFactors=FALSE)
EnsembIDs2 <- read.table(file="/public/ZhangJinLab/backdata/lahdata/M2statgenes.txt",sep="\t",header=FALSE,stringsAsFactors=FALSE)
EnsembIDs <- rbind(EnsembIDs1,EnsembIDs2)

mapInfo <- read.table(file="/public/ZhangJinLab/backdata/lahdata/AA/features.tsv.gz",sep="\t",header=FALSE,stringsAsFactors=FALSE)
matchIndexes <- match(EnsembIDs[,2],mapInfo[,1])
geneVec <- mapInfo[matchIndexes[which(!is.na(matchIndexes))],2]

res <- c()
geneIDs <- c()
for (genename in geneVec){
	if(length(which(rownames(normDataAA)==genename)) == 1){
		AA <- normDataAA[which(rownames(normDataAA)==genename),]
		AA <- mean(AA)

		BB <- normDataBB[which(rownames(normDataBB)==genename),]
		BB <- mean(BB)

		DD <- normDataDD[which(rownames(normDataDD)==genename),]
		DD <- mean(DD)

		if(length(res) == 0){
			res <- c(AA,BB,DD)
			geneIDs <- c(geneIDs,genename)
		}else{
			res <- rbind(res,c(AA,BB,DD))
			geneIDs <- c(geneIDs,genename)
		}
	}
}
rownames(res) <- geneIDs
colnames(res) <- c("AA","BB","DD")
res <- res[which(rowSums(res) > 0),]

library(pheatmap)
pdf(file="/public/ZhangJinLab/backdata/lahdata/M1M2statgenesheatmap.pdf",width=5,height=10)
pheatmap(res,scale="row",cluster_rows=TRUE,cluster_cols=FALSE,show_rownames=TRUE,show_colnames=TRUE,color = colorRampPalette(c("#6B70B0","white","red"))(50),fontsize=20)
dev.off()


res <- c()
geneIDs <- c()
for (genename in geneVec){
	if(length(which(rownames(normDataAA)==genename)) == 1){
		AA <- normDataAA[which(rownames(normDataAA)==genename),]
		CD80AA <- normDataAA["CD80",]
		AACD80 <- mean(AA[which(CD80AA > 0)])
		AANonCD80 <- mean(AA[which(CD80AA == 0)])	

		BB <- normDataBB[which(rownames(normDataBB)==genename),]
		CD80BB <- normDataBB["CD80",]
		BBCD80 <- mean(BB[which(CD80BB > 0)])
		BBNonCD80 <- mean(BB[which(CD80BB == 0)])

		DD <- normDataDD[which(rownames(normDataDD)==genename),]
		CD80DD <- normDataDD["CD80",]
		DDCD80 <- mean(DD[which(CD80DD > 0)])
		DDNonCD80 <- mean(DD[which(CD80DD == 0)])
		
		if(length(res) == 0){
			res <- c(AACD80,BBCD80,DDCD80)
			geneIDs <- c(geneIDs,genename)
		}else{
			res <- rbind(res,c(AACD80,BBCD80,DDCD80))
			geneIDs <- c(geneIDs,genename)
		}
	}
}

rownames(res) <- geneIDs
colnames(res) <- c("AACD80","BBCD80","DDCD80")
res <- res[which(rowSums(res) > 0),]

library(pheatmap)
pdf(file="/public/ZhangJinLab/backdata/lahdata/M1M2statgenesCD80heatmap.pdf",width=5,height=10)
pheatmap(res,scale="row",show_rownames=TRUE,show_colnames=TRUE,cluster_rows=TRUE,cluster_cols=FALSE,color = colorRampPalette(c("#6B70B0","white","red"))(100),fontsize=20)
dev.off()
