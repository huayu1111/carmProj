library(DESeq2)
library(ggrepel)

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

expDataAA_r1 <- rowSums(rawDataAA[,sample(ncol(rawDataAA),10000)])
expDataAA_r2 <- rowSums(rawDataAA[,sample(ncol(rawDataAA),10000)])
expDataBB_r1 <- rowSums(rawDataBB[,sample(ncol(rawDataBB),10000)])
expDataBB_r2 <- rowSums(rawDataBB[,sample(ncol(rawDataBB),10000)])
expDataDD_r1 <- rowSums(rawDataDD[,sample(ncol(rawDataDD),10000)])
expDataDD_r2 <- rowSums(rawDataDD[,sample(ncol(rawDataDD),10000)])

matchIndexes <- match(rownames(rawDataAA),rownames(rawDataBB))
expDataAABB <- cbind(expDataAA_r1[which(!is.na(matchIndexes))],expDataAA_r2[which(!is.na(matchIndexes))],expDataBB_r1[matchIndexes[which(!is.na(matchIndexes))]],expDataBB_r2[matchIndexes[which(!is.na(matchIndexes))]])

matchIndexes <- match(rownames(rawDataAA),rownames(rawDataDD))
expDataAADD <- cbind(expDataAA_r1[which(!is.na(matchIndexes))],expDataAA_r2[which(!is.na(matchIndexes))],expDataDD_r1[matchIndexes[which(!is.na(matchIndexes))]],expDataDD_r2[matchIndexes[which(!is.na(matchIndexes))]])

matchIndexes <- match(rownames(rawDataBB),rownames(rawDataDD))
expDataBBDD <- cbind(expDataBB_r1[which(!is.na(matchIndexes))],expDataBB_r2[which(!is.na(matchIndexes))],expDataDD_r1[matchIndexes[which(!is.na(matchIndexes))]],expDataDD_r2[matchIndexes[which(!is.na(matchIndexes))]])

colnames(expDataAABB[,c(1,2,3,4)]) <- c("AAR1","AAR2","BBR1","BBR2")
conditions <- factor(c(rep("AA",2),rep("BB",2)),levels=c("AA","BB"))
coldata <- data.frame(row.names = colnames(expDataAABB[,c(1,2,3,4)]), conditions)
dds <- DESeqDataSetFromMatrix(countData=expDataAABB[,c(1,2,3,4)], colData=coldata, design=~conditions)
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj),]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
matchIndexes <- match(resdata[,1],mapData[,1])
rownames(resdata) <- resdata[,1]
resdata[,1] <- mapData[matchIndexes,2]
colnames(resdata)[1] <- c("gene_name")

targetedGenes <- c("IL18","CCL15","CD7","CD48","TNFAIP2","IL15","IL1A","TNFSF13B","CCL8","IL1B","IL23A","CXCL2","CXCL10","IL32","CCL19","CXCL11","CXCL9","CD83","CD82","IL1RN","IL2RA","TLR8","TNFSF10","IL27","CD80","CD40","IL4I1","CD38","IL15RA","CD36","CD93","CSF1R","TNFRSF11A","TREM2","NPC1","CD209","CCL13","CD200R1","IL10","CD180","LPL","CD1E","CCL26","ARG1")
nonnaIndexes <- which(!is.na(resdata$pvalue))
resdata <- resdata[nonnaIndexes,]
matchIndexes <- match(targetedGenes,resdata$gene_name)
resdata$classlabels <- rep("Nonsignificance",nrow(resdata))
resdata$classlabels[which(resdata$pvalue < 0.05 & resdata$log2FoldChange > 0.5)] <- "Upregulation"
resdata$classlabels[which(resdata$pvalue < 0.05 & resdata$log2FoldChange < -0.5)] <- "Downregulation"
resdata$label <-''
resdata$label[matchIndexes] <- targetedGenes

resdata$size <- rep(0.1,nrow(resdata))
resdata$pvalue <- -log10(resdata$pvalue)
infIndexes <- which(is.infinite(resdata$pvalue))
noninfIndexes <- which(!is.infinite(resdata$pvalue))
maxValue <- max(resdata$pvalue[noninfIndexes])
resdata$pvalue[infIndexes] <- maxValue

labeldata <- resdata[matchIndexes,]
labeldata <- labeldata[which(labeldata$pvalue > -log10(0.05) & abs(labeldata$log2FoldChange) > 0.5),]
pdf("/home/Yuhua/lahmerge/cdata//AA_BB_Diff_Markers_Volcano.pdf",width=5,height=5)
p <- ggplot(resdata,aes(x=log2FoldChange,y=pvalue,colour=factor(resdata$classlabels,levels=c("Downregulation","Nonsignificance","Upregulation")))) + geom_point(size=resdata$size,shape=20,show.legend = F) + scale_colour_manual(values=c("blue","gray","red")) + geom_text_repel(data=labeldata, aes(x=log2FoldChange,y=pvalue,label=label), colour="black", size=2,direction="both",min.segment.length = 0.05,segment.alpha=0.6,max.overlaps =300,nudge_x = 0.2,nudge_y=0.2)
p <- p + xlab(expression(log[2]("fold change"))) + ylab(expression(-log[10]("p-value")))
p <- p + theme(legend.title=element_blank())
p <- p + scale_y_continuous(limits = c(0,330))
p <- p + scale_x_continuous(limits = c(-4,4))
p <- p + annotate(geom="text", x=-3.5, y=250, label=paste("Down:",length(which(resdata$classlabels=="Downregulation")),sep=""), color="blue")
p <- p + annotate(geom="text", x=3.5, y=250, label=paste("Up:",length(which(resdata$classlabels=="Upregulation")),sep=""), color="red")
p <- p + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11))
print(p)
dev.off()


colnames(expDataAADD[,c(1,2,3,4)]) <- c("AAR1","AAR2","DDR1","DDR2")
conditions <- factor(c(rep("AA",2),rep("DD",2)),levels=c("AA","DD"))
coldata <- data.frame(row.names = colnames(expDataAADD[,c(1,2,3,4)]), conditions)
dds <- DESeqDataSetFromMatrix(countData=expDataAADD[,c(1,2,3,4)], colData=coldata, design=~conditions)
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj),]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
matchIndexes <- match(resdata[,1],mapData[,1])
rownames(resdata) <- resdata[,1]
resdata[,1] <- mapData[matchIndexes,2]
colnames(resdata)[1] <- c("gene_name")

targetedGenes <- c("IL18","CCL15","CD7","CD48","TNFAIP2","IL15","IL1A","TNFSF13B","CCL8","IL1B","IL23A","CXCL2","CXCL10","IL32","CCL19","CXCL11","CXCL9","CD83","CD82","IL1RN","IL2RA","TLR8","TNFSF10","IL27","CD80","CD40","IL4I1","CD38","IL15RA","CD36","CD93","CSF1R","TNFRSF11A","TREM2","NPC1","CD209","CCL13","CD200R1","IL10","CD180","LPL","CD1E","CCL26","ARG1")
nonnaIndexes <- which(!is.na(resdata$pvalue))
resdata <- resdata[nonnaIndexes,]
matchIndexes <- match(targetedGenes,resdata$gene_name)
resdata$classlabels <- rep("Nonsignificance",nrow(resdata))
resdata$classlabels[which(resdata$pvalue < 0.05 & resdata$log2FoldChange > 0.5)] <- "Upregulation"
resdata$classlabels[which(resdata$pvalue < 0.05 & resdata$log2FoldChange < -0.5)] <- "Downregulation"
resdata$label <-''
resdata$label[matchIndexes] <- targetedGenes

resdata$size <- rep(0.1,nrow(resdata))
resdata$pvalue <- -log10(resdata$pvalue)
infIndexes <- which(is.infinite(resdata$pvalue))
noninfIndexes <- which(!is.infinite(resdata$pvalue))
maxValue <- max(resdata$pvalue[noninfIndexes])
resdata$pvalue[infIndexes] <- maxValue

labeldata <- resdata[matchIndexes,]
labeldata <- labeldata[which(labeldata$pvalue > -log10(0.05) & abs(labeldata$log2FoldChange) > 0.5),]
pdf("/home/Yuhua/lahmerge/cdata//AA_DD_Diff_Markers_Volcano.pdf",width=5,height=5)
p <- ggplot(resdata,aes(x=log2FoldChange,y=pvalue,colour=factor(resdata$classlabels,levels=c("Downregulation","Nonsignificance","Upregulation")))) + geom_point(size=resdata$size,shape=20,show.legend = F) + scale_colour_manual(values=c("blue","gray","red")) + geom_text_repel(data=labeldata, aes(x=log2FoldChange,y=pvalue,label=label), colour="black", size=2,direction="both",min.segment.length = 0.05,segment.alpha=0.6,max.overlaps =300,nudge_x = 0.2,nudge_y=0.2)
p <- p + xlab(expression(log[2]("fold change"))) + ylab(expression(-log[10]("p-value")))
p <- p + theme(legend.title=element_blank())
p <- p + scale_y_continuous(limits = c(0,330))
p <- p + scale_x_continuous(limits = c(-4,4))
p <- p + annotate(geom="text", x=-3.5, y=250, label=paste("Down:",length(which(resdata$classlabels=="Downregulation")),sep=""), color="blue")
p <- p + annotate(geom="text", x=3.5, y=250, label=paste("Up:",length(which(resdata$classlabels=="Upregulation")),sep=""), color="red")
p <- p + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11))
print(p)
dev.off()


colnames(expDataBBDD[,c(1,2,3,4)]) <- c("BBR1","BBR2","DDR1","DDR2")
conditions <- factor(c(rep("BB",2),rep("DD",2)),levels=c("BB","DD"))
coldata <- data.frame(row.names = colnames(expDataBBDD[,c(1,2,3,4)]), conditions)
dds <- DESeqDataSetFromMatrix(countData=expDataBBDD[,c(1,2,3,4)], colData=coldata, design=~conditions)
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj),]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
matchIndexes <- match(resdata[,1],mapData[,1])
rownames(resdata) <- resdata[,1]
resdata[,1] <- mapData[matchIndexes,2]
colnames(resdata)[1] <- c("gene_name")

targetedGenes <- c("IL18","CCL15","CD7","CD48","TNFAIP2","IL15","IL1A","TNFSF13B","CCL8","IL1B","IL23A","CXCL2","CXCL10","IL32","CCL19","CXCL11","CXCL9","CD83","CD82","IL1RN","IL2RA","TLR8","TNFSF10","IL27","CD80","CD40","IL4I1","CD38","IL15RA","CD36","CD93","CSF1R","TNFRSF11A","TREM2","NPC1","CD209","CCL13","CD200R1","IL10","CD180","LPL","CD1E","CCL26","ARG1")
nonnaIndexes <- which(!is.na(resdata$pvalue))
resdata <- resdata[nonnaIndexes,]
matchIndexes <- match(targetedGenes,resdata$gene_name)
resdata$classlabels <- rep("Nonsignificance",nrow(resdata))
resdata$classlabels[which(resdata$pvalue < 0.05 & resdata$log2FoldChange > 0.5)] <- "Upregulation"
resdata$classlabels[which(resdata$pvalue < 0.05 & resdata$log2FoldChange < -0.5)] <- "Downregulation"
resdata$label <-''
resdata$label[matchIndexes] <- targetedGenes

resdata$size <- rep(0.1,nrow(resdata))
resdata$pvalue <- -log10(resdata$pvalue)
infIndexes <- which(is.infinite(resdata$pvalue))
noninfIndexes <- which(!is.infinite(resdata$pvalue))
maxValue <- max(resdata$pvalue[noninfIndexes])
resdata$pvalue[infIndexes] <- maxValue

labeldata <- resdata[matchIndexes,]
labeldata <- labeldata[which(labeldata$pvalue > -log10(0.05) & abs(labeldata$log2FoldChange) > 0.5),]
pdf("/home/Yuhua/lahmerge/cdata//BB_DD_Diff_Markers_Volcano.pdf",width=5,height=5)
p <- ggplot(resdata,aes(x=log2FoldChange,y=pvalue,colour=factor(resdata$classlabels,levels=c("Downregulation","Nonsignificance","Upregulation")))) + geom_point(size=resdata$size,shape=20,show.legend = F) + scale_colour_manual(values=c("blue","gray","red")) + geom_text_repel(data=labeldata, aes(x=log2FoldChange,y=pvalue,label=label), colour="black", size=2,direction="both",min.segment.length = 0.05,segment.alpha=0.6,max.overlaps =300,nudge_x = 0.2,nudge_y=0.2)
p <- p + xlab(expression(log[2]("fold change"))) + ylab(expression(-log[10]("p-value")))
p <- p + theme(legend.title=element_blank())
p <- p + scale_y_continuous(limits = c(0,330))
p <- p + scale_x_continuous(limits = c(-4,4))
p <- p + annotate(geom="text", x=-3.5, y=250, label=paste("Down:",length(which(resdata$classlabels=="Downregulation")),sep=""), color="blue")
p <- p + annotate(geom="text", x=3.5, y=250, label=paste("Up:",length(which(resdata$classlabels=="Upregulation")),sep=""), color="red")
p <- p + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11))
print(p)
dev.off()