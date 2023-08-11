tpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

mapData <- read.table(file="/home/Yuhua/lahmerge/cdata/features.tsv.gz",sep="\t")
rawDataAA <- read.table(file="/home/Yuhua/lahmerge/cdata/AA/10X_LAH_InHouse_rawcounts.txt",row.names=1,header=TRUE)
matchIndexes <- match(rownames(rawDataAA),mapData[,2])
rawDataAA <- rawDataAA[!is.na(matchIndexes),]
rownames(rawDataAA) <- mapData[matchIndexes[!is.na(matchIndexes)],1]
rawDataBB <- read.table(file="/home/Yuhua/lahmerge/cdata/BB/10X_LAH_InHouse_rawcounts.txt",row.names=1,header=TRUE)
matchIndexes <- match(rownames(rawDataBB),mapData[,2])
rawDataBB <- rawDataBB[!is.na(matchIndexes),]
rownames(rawDataBB) <- mapData[matchIndexes[!is.na(matchIndexes)],1]
rawDataDD <- read.table(file="/home/Yuhua/lahmerge/cdata/DD/10X_LAH_InHouse_rawcounts.txt",row.names=1,header=TRUE)
matchIndexes <- match(rownames(rawDataDD),mapData[,2])
rawDataDD <- rawDataDD[!is.na(matchIndexes),]
rownames(rawDataDD) <- mapData[matchIndexes[!is.na(matchIndexes)],1]

expDataAA <- rowSums(rawDataAA)
expDataBB <- rowSums(rawDataBB)
expDataDD <- rowSums(rawDataDD)
matchIndexes <- match(rownames(rawDataAA),rownames(rawDataBB))
expDataAABB <- cbind(expDataAA[which(!is.na(matchIndexes))],expDataBB[matchIndexes[which(!is.na(matchIndexes))]])
matchIndexes <- match(rownames(expDataAABB),rownames(rawDataDD))
expData <- cbind(expDataAABB[which(!is.na(matchIndexes)),],expDataDD[matchIndexes[which(!is.na(matchIndexes))]])

colnames(expData) <- c("AA","BB","DD")
genemapData <- read.table(file="/home/Yuhua/lahmerge/cdata/gencode.v28.genetype.tab",sep="\t",stringsAsFactors=FALSE,header=FALSE)
genelenData <- read.table(file="/home/Yuhua/lahmerge/cdata/gene_cds_length.txt",sep=" ",stringsAsFactors=FALSE,header=FALSE)
matchIndexes <- match(rownames(expData),gsub("\\.\\d+","",genelenData[,1]))
expData <- expData[which(!is.na(matchIndexes)),]

expData <- tpm(expData,as.numeric(genelenData[matchIndexes[which(!is.na(matchIndexes))],2]))
matchIndexes <- match(rownames(expData),gsub("\\.\\d+","",genemapData[,1]))
expData <- cbind(genemapData[matchIndexes,3],expData)
write.table(expData,file="/home/Yuhua/lahmerge/cdata/AABBDD_gene_TPM.txt",sep='\t',quote=FALSE,row.names=T,col.names=T)


expDataAACD80 <- rowSums(rawDataAA[,which(rawDataAA["ENSG00000121594",] > 0)])
expDataAANonCD80 <- rowSums(rawDataAA[,which(rawDataAA["ENSG00000121594",] == 0)])
expDataBBCD80 <- rowSums(rawDataBB[,which(rawDataBB["ENSG00000121594",] > 0)])
expDataBBNonCD80 <- rowSums(rawDataBB[,which(rawDataBB["ENSG00000121594",] == 0)])
expDataDDCD80 <- rowSums(rawDataDD[,which(rawDataDD["ENSG00000121594",] > 0)])
expDataDDNonCD80 <- rowSums(rawDataDD[,which(rawDataDD["ENSG00000121594",] == 0)])
matchIndexes <- match(rownames(rawDataAA),rownames(rawDataBB))
expDataAABB <- cbind(expDataAACD80[which(!is.na(matchIndexes))],expDataAANonCD80[which(!is.na(matchIndexes))],expDataBBCD80[matchIndexes[which(!is.na(matchIndexes))]],expDataBBNonCD80[matchIndexes[which(!is.na(matchIndexes))]])
matchIndexes <- match(rownames(expDataAABB),rownames(rawDataDD))
expData <- cbind(expDataAABB[which(!is.na(matchIndexes)),],expDataDDCD80[matchIndexes[which(!is.na(matchIndexes))]],expDataDDNonCD80[matchIndexes[which(!is.na(matchIndexes))]])

colnames(expData) <- c("AACD80","AANonCD80","BBCD80","BBNonCD80","DDCD80","DDNonCD80")
genemapData <- read.table(file="/home/Yuhua/lahmerge/cdata/gencode.v28.genetype.tab",sep="\t",stringsAsFactors=FALSE,header=FALSE)
genelenData <- read.table(file="/home/Yuhua/lahmerge/cdata/gene_cds_length.txt",sep=" ",stringsAsFactors=FALSE,header=FALSE)
matchIndexes <- match(rownames(expData),gsub("\\.\\d+","",genelenData[,1]))
expData <- expData[which(!is.na(matchIndexes)),]

expData <- tpm(expData,as.numeric(genelenData[matchIndexes[which(!is.na(matchIndexes))],2]))
matchIndexes <- match(rownames(expData),gsub("\\.\\d+","",genemapData[,1]))
expData <- cbind(genemapData[matchIndexes,3],expData)
write.table(expData,file="/home/Yuhua/lahmerge/cdata/AABBDDCD80_gene_TPM.txt",sep='\t',quote=FALSE,row.names=T,col.names=T)
