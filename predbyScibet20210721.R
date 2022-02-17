library(scibet)
library(tidyverse)
library(ggplot2)
library(viridis)
library(Matrix)

tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

genemapData <- read.table(file="/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/iNK/gencode.v28.genetype.tab",sep="\t",stringsAsFactors=FALSE,header=FALSE)
genelenData <- read.table(file="/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/iNK/gene_cds_length.txt",sep=" ",stringsAsFactors=FALSE,header=FALSE)
matchIndexes <- match(genelenData[,1],genemapData[,1])
lenInfo <- cbind(genemapData[matchIndexes[which(!is.na(matchIndexes))],3],genelenData[which(!is.na(matchIndexes)),2])

files <- Sys.glob(file.path("/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/iNK/pubdata/","*","*","GRCh38","matrix.mtx"))

pubdataset <- c()
for (filename in files){
	mydata_total <- readMM(file = filename)
	featurename <- gsub("matrix.mtx","genes.tsv",filename)
	cellname <- gsub("matrix.mtx","barcodes.tsv",filename)
	mydata.feature.names = read.delim(featurename,
                           header = FALSE,
                           stringsAsFactors = FALSE)
	mydata.barcode.names = read.delim(cellname,
                           header = FALSE,
                           stringsAsFactors = FALSE)
	colnames(mydata_total) = gsub("\\-1","",mydata.barcode.names$V1)
	rownames(mydata_total) = mydata.feature.names$V2
	posIndex <- regexpr('pubdata\\/..*?\\/',filename)
	sampleid <- substring(filename,posIndex+9,posIndex+attr(posIndex,'match.length')-2)
	cat(as.character(sampleid),"\n")
	mapfile <- paste("/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/iNK/pubdata/",sampleid,"/",sampleid,"/",sampleid,".csv",sep="")
	mapdata <- read.csv(mapfile,header=TRUE,stringsAsFactors=FALSE)
	mapdata <- mapdata[!is.na(mapdata[,2]),]
	matchIndexes <- match(colnames(mydata_total),mapdata[,1])
	mydata_total <- mydata_total[,which(!is.na(matchIndexes))]
	labelnames <- mapdata[matchIndexes[which(!is.na(matchIndexes))],2]
	matchIndexes <- match(rownames(mydata_total),lenInfo[,1])
	mydata_total <- mydata_total[which(!is.na(matchIndexes)),]
	lenVec <- lenInfo[matchIndexes[which(!is.na(matchIndexes))],2]
	mydata_total <- tpm3(mydata_total,as.numeric(lenVec))
	
	mydata_total <- t(mydata_total)
	mydata_total <- cbind(as.data.frame(as.matrix(mydata_total)),as.data.frame(labelnames))
	colnames(mydata_total)[ncol(mydata_total)] <- "label"
	rownames(mydata_total) <- NULL
	
	if(length(pubdataset)==0){
		pubdataset <- mydata_total
	}else{
		pubdataset <- 	rbind(pubdataset,mydata_total)
	}
}

rawData <- read.table(file="/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/AA/10X_LAH_InHouse_rawcounts.txt",row.names=1,header=TRUE)
normData <- read.table(file="/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/AA/10X_LAH_InHouse_normcounts.txt",row.names=1,header=TRUE)
cellmatchIndexes <- match(colnames(normData),colnames(rawData))
rawData <- rawData[,cellmatchIndexes]

matchIndexes <- match(rownames(rawData),lenInfo[,1])
countData <- rawData[which(!is.na(matchIndexes)),]
lenVec <- lenInfo[matchIndexes[which(!is.na(matchIndexes))],2]
tpmData <- tpm3(countData,as.numeric(lenVec))

matchIndexes <- match(rownames(tpmData),colnames(pubdataset))
tpmData <- tpmData[which(!is.na(matchIndexes)),]

pubdataset <- pubdataset[,c(matchIndexes[which(!is.na(matchIndexes))],ncol(pubdataset))]

orgin <- read.table(file="/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/AA/10X_LAH_InHouse_cellidentity.txt",header=FALSE,stringsAsFactors=FALSE)[,1]
orgin <- orgin + 1
orgin <- paste("C",orgin,sep="")

prd <- SciBet(pubdataset, t(tpmData))
save(orgin,prd,file="/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/AA/predinfo.RData")
pdf("/public/ZhangJinLab/project_erv/Dpan/IPSAnaly/AA/predcelltypes_heatmap.pdf",width=10,height=10)
Confusion_heatmap(orgin, prd)
dev.off()
