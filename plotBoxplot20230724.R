library(Seurat)
library(ggplot2)
library(beeswarm)

load("/public/ZhangJinLab/backdata/lahdata/AA/10X_LAH_InHouse.RData")
normDataAA <- GetAssayData(object = mydata_total[["RNA"]], slot = "data")
load("/public/ZhangJinLab/backdata/lahdata/BB/10X_LAH_InHouse.RData")
normDataBB <- GetAssayData(object = mydata_total[["RNA"]], slot = "data")
load("/public/ZhangJinLab/backdata/lahdata/DD/10X_LAH_InHouse.RData")
normDataDD <- GetAssayData(object = mydata_total[["RNA"]], slot = "data")

# Part 1
EnsembIDs <- read.table(file="/public/ZhangJinLab/backdata/lahdata/M2genes",sep="\t",header=FALSE,stringsAsFactors=FALSE)
mapInfo <- read.table(file="/public/ZhangJinLab/backdata/lahdata/AA/features.tsv.gz",sep="\t",header=FALSE,stringsAsFactors=FALSE)
matchIndexes <- match(EnsembIDs[,2],mapInfo[,1])
geneVec <- mapInfo[matchIndexes[which(!is.na(matchIndexes))],2]
EnsembIDs <- EnsembIDs[which(!is.na(matchIndexes)),1]

res <- c()
iterer=1
for (genename in geneVec){
	if(length(which(rownames(normDataAA)==genename)) == 1){
		AA <- normDataAA[which(rownames(normDataAA)==genename),]
		AA <- AA[which(AA > 0)]
		cat(length(AA)/ncol(normDataAA),"\n")

		BB <- normDataBB[which(rownames(normDataBB)==genename),]
		BB <- BB[which(BB > 0)]
		cat(length(BB)/ncol(normDataBB),"\n")

		DD <- normDataDD[which(rownames(normDataDD)==genename),]
		DD <- DD[which(DD > 0)]
		cat(length(DD)/ncol(normDataDD),"\n")

		a <- wilcox.test(AA,BB)
		b <- wilcox.test(AA,DD)
		c <- wilcox.test(BB,DD)
		if(length(res) == 0){
			res <- c(genename,(length(AA)/ncol(normDataAA)),(length(BB)/ncol(normDataBB)),(length(DD)/ncol(normDataDD)),a$p.value,b$p.value,c$p.value)
		}else{
			res <- rbind(res,c(genename,(length(AA)/ncol(normDataAA)),(length(BB)/ncol(normDataBB)),(length(DD)/ncol(normDataDD)),a$p.value,b$p.value,c$p.value))
		}

		exp_des_label <- rep(c("AA","BB","DD"),c(length(AA),length(BB),length(DD)))
		exp_des_value <- as.numeric(c(AA,BB,DD))
		exp_des <- data.frame(cbind(exp_des_label,exp_des_value))
		colnames(exp_des) <- c("Label","Value")
		exp_des$Label=factor(exp_des_label,levels=c("AA","BB","DD"))
		exp_des$value <- as.numeric(exp_des_value)

		library(colorspace)
		q4 <- qualitative_hcl(4, palette = "Pastel 1")
		library(tidyverse)
		library(ggplot2)
		library(ggbeeswarm)
		df_means <- exp_des %>% group_by(Label) %>% summarise(value=mean(value))

		pdf(paste("/public/ZhangJinLab/backdata/lahdata/M2Gene/",genename,"_boxplot.pdf",sep=""),width=8,height=10)
		p <- ggplot(data=exp_des,aes(x=Label,y=value)) + geom_quasirandom(aes(color=Label),method="frowney",width=0.3,size=1) + scale_color_manual(values=q4[c(1,2,4)]) + geom_boxplot(width=0.05,outlier.shape=NA,color="grey",fill="white",linetype=2) + stat_summary(fun.y=mean, geom="point", shape=20, size=10,  alpha=0.7, color="#8B814C")
		p <- p + xlab("") + ylab("Normalized expression level") + ggtitle("") + theme_classic()
		p <- p + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),panel.background=element_blank(), rect=element_rect(linetype=1), axis.line=element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11),legend.position="none") + theme(axis.text.x=element_text(vjust = 0.6, angle = 45))
		print(p)
		dev.off()
		iterer <- iterer + 1
	}
}
write.table(res,file="/public/ZhangJinLab/backdata/lahdata/M2Gene/M2genes_stat.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

res <- c()
iterer=1
for (genename in geneVec){
	if(length(which(rownames(normDataAA)==genename)) == 1){
		AA <- normDataAA[which(rownames(normDataAA)==genename),]
		CD80AA <- normDataAA["CD80",]
		AA <- AA[which(CD80AA > 0)]
		AA <- AA[which(AA > 0)]
		
		cat(length(AA)/ncol(normDataAA),"\n")

		BB <- normDataBB[which(rownames(normDataBB)==genename),]
		CD80BB <- normDataBB["CD80",]
		BB <- BB[which(CD80BB > 0)]
		BB <- BB[which(BB > 0)]
		cat(length(BB)/ncol(normDataBB),"\n")

		DD <- normDataDD[which(rownames(normDataDD)==genename),]
		CD80DD <- normDataDD["CD80",]
		DD <- DD[which(CD80DD > 0)]
		DD <- DD[which(DD > 0)]
		cat(length(DD)/ncol(normDataDD),"\n")
		
		if(length(AA) > 3 && length(BB) > 3 && length(DD) > 3){
			a <- wilcox.test(AA,BB)
			b <- wilcox.test(AA,DD)
			c <- wilcox.test(BB,DD)
			if(length(res) == 0){
				res <- c(genename,(length(AA)/ncol(normDataAA)),(length(BB)/ncol(normDataBB)),(length(DD)/ncol(normDataDD)),a$p.value,b$p.value,c$p.value)
			}else{
				res <- rbind(res,c(genename,(length(AA)/ncol(normDataAA)),(length(BB)/ncol(normDataBB)),(length(DD)/ncol(normDataDD)),a$p.value,b$p.value,c$p.value))
			}

			exp_des_label <- rep(c("AA","BB","DD"),c(length(AA),length(BB),length(DD)))
			exp_des_value <- as.numeric(c(AA,BB,DD))
			exp_des <- data.frame(cbind(exp_des_label,exp_des_value))
			colnames(exp_des) <- c("Label","Value")
			exp_des$Label=factor(exp_des_label,levels=c("AA","BB","DD"))
			exp_des$value <- as.numeric(exp_des_value)

			library(colorspace)
			q4 <- qualitative_hcl(4, palette = "Pastel 1")
			library(tidyverse)
			library(ggplot2)
			library(ggbeeswarm)
			df_means <- exp_des %>% group_by(Label) %>% summarise(value=mean(value))

			pdf(paste("/public/ZhangJinLab/backdata/lahdata/M2Gene/",genename,"_CD80_boxplot.pdf",sep=""),width=8,height=10)
			p <- ggplot(data=exp_des, aes(x=Label,y=value)) + geom_quasirandom(aes(color=Label),method="frowney",width=0.3,size=1) + scale_color_manual(values=q4[c(1,2,4)]) + geom_boxplot(width=0.05,outlier.shape=NA,color="grey",fill="white",linetype=2) + stat_summary(fun.y=mean, geom="point", shape=20, size=10,  alpha=0.7, color="#8B814C")
			p <- p + xlab("") + ylab("Normalized expression level") + ggtitle("") + theme_classic()
			p <- p + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),panel.background=element_blank(), rect=element_rect(linetype=1), axis.line=element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11),legend.position="none") + theme(axis.text.x=element_text(vjust = 0.6, angle = 45))
			print(p)
			dev.off()
			iterer <- iterer + 1
		}
	}
}
write.table(res,file="/public/ZhangJinLab/backdata/lahdata/M2Gene/M2genes_CD80_stat.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)


# Part 2
EnsembIDs <- read.table(file="/public/ZhangJinLab/backdata/lahdata/Tgenes",sep="\t",header=FALSE,stringsAsFactors=FALSE)
mapInfo <- read.table(file="/public/ZhangJinLab/backdata/lahdata/AA/features.tsv.gz",sep="\t",header=FALSE,stringsAsFactors=FALSE)
matchIndexes <- match(EnsembIDs[,2],mapInfo[,1])
geneVec <- mapInfo[matchIndexes[which(!is.na(matchIndexes))],2]
EnsembIDs <- EnsembIDs[which(!is.na(matchIndexes)),1]

res <- c()
iterer=1
for (genename in geneVec){
	if(length(which(rownames(normDataAA)==genename)) == 1){
		AA <- normDataAA[which(rownames(normDataAA)==genename),]
		AA <- AA[which(AA > 0)]
		cat(length(AA)/ncol(normDataAA),"\n")

		BB <- normDataBB[which(rownames(normDataBB)==genename),]
		BB <- BB[which(BB > 0)]
		cat(length(BB)/ncol(normDataBB),"\n")

		DD <- normDataDD[which(rownames(normDataDD)==genename),]
		DD <- DD[which(DD > 0)]
		cat(length(DD)/ncol(normDataDD),"\n")

		a <- wilcox.test(AA,BB)
		b <- wilcox.test(AA,DD)
		c <- wilcox.test(BB,DD)
		if(length(res) == 0){
			res <- c(genename,(length(AA)/ncol(normDataAA)),(length(BB)/ncol(normDataBB)),(length(DD)/ncol(normDataDD)),a$p.value,b$p.value,c$p.value)
		}else{
			res <- rbind(res,c(genename,(length(AA)/ncol(normDataAA)),(length(BB)/ncol(normDataBB)),(length(DD)/ncol(normDataDD)),a$p.value,b$p.value,c$p.value))
		}

		exp_des_label <- rep(c("AA","BB","DD"),c(length(AA),length(BB),length(DD)))
		exp_des_value <- as.numeric(c(AA,BB,DD))
		exp_des <- data.frame(cbind(exp_des_label,exp_des_value))
		colnames(exp_des) <- c("Label","Value")
		exp_des$Label=factor(exp_des_label,levels=c("AA","BB","DD"))
		exp_des$value <- as.numeric(exp_des_value)

		library(colorspace)
		q4 <- qualitative_hcl(4, palette = "Pastel 1")
		library(tidyverse)
		library(ggplot2)
		library(ggbeeswarm)
		df_means <- exp_des %>% group_by(Label) %>% summarise(value=mean(value))

		pdf(paste("/public/ZhangJinLab/backdata/lahdata/TGene/",genename,"_boxplot.pdf",sep=""),width=8,height=10)
		p <- ggplot(data=exp_des,aes(x=Label,y=value)) + geom_quasirandom(aes(color=Label),method="frowney",width=0.3,size=1) + scale_color_manual(values=q4[c(1,2,4)]) + geom_boxplot(width=0.05,outlier.shape=NA,color="grey",fill="white",linetype=2) + stat_summary(fun.y=mean, geom="point", shape=20, size=10,  alpha=0.7, color="#8B814C")
		p <- p + xlab("") + ylab("Normalized expression level") + ggtitle("") + theme_classic()
		p <- p + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),panel.background=element_blank(), rect=element_rect(linetype=1), axis.line=element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11),legend.position="none") + theme(axis.text.x=element_text(vjust = 0.6, angle = 45))
		print(p)
		dev.off()
		iterer <- iterer + 1
	}
}
write.table(res,file="/public/ZhangJinLab/backdata/lahdata/TGene/Tgenes_stat.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

res <- c()
iterer=1
for (genename in geneVec){
	if(length(which(rownames(normDataAA)==genename)) == 1){
		AA <- normDataAA[which(rownames(normDataAA)==genename),]
		CD80AA <- normDataAA["CD80",]
		AA <- AA[which(CD80AA > 0)]
		AA <- AA[which(AA > 0)]
		
		cat(length(AA)/ncol(normDataAA),"\n")

		BB <- normDataBB[which(rownames(normDataBB)==genename),]
		CD80BB <- normDataBB["CD80",]
		BB <- BB[which(CD80BB > 0)]
		BB <- BB[which(BB > 0)]
		cat(length(BB)/ncol(normDataBB),"\n")

		DD <- normDataDD[which(rownames(normDataDD)==genename),]
		CD80DD <- normDataDD["CD80",]
		DD <- DD[which(CD80DD > 0)]
		DD <- DD[which(DD > 0)]
		cat(length(DD)/ncol(normDataDD),"\n")
		
		if(length(AA) > 3 && length(BB) > 3 && length(DD) > 3){
			a <- wilcox.test(AA,BB)
			b <- wilcox.test(AA,DD)
			c <- wilcox.test(BB,DD)
			if(length(res) == 0){
				res <- c(genename,(length(AA)/ncol(normDataAA)),(length(BB)/ncol(normDataBB)),(length(DD)/ncol(normDataDD)),a$p.value,b$p.value,c$p.value)
			}else{
				res <- rbind(res,c(genename,(length(AA)/ncol(normDataAA)),(length(BB)/ncol(normDataBB)),(length(DD)/ncol(normDataDD)),a$p.value,b$p.value,c$p.value))
			}

			exp_des_label <- rep(c("AA","BB","DD"),c(length(AA),length(BB),length(DD)))
			exp_des_value <- as.numeric(c(AA,BB,DD))
			exp_des <- data.frame(cbind(exp_des_label,exp_des_value))
			colnames(exp_des) <- c("Label","Value")
			exp_des$Label=factor(exp_des_label,levels=c("AA","BB","DD"))
			exp_des$value <- as.numeric(exp_des_value)

			library(colorspace)
			q4 <- qualitative_hcl(4, palette = "Pastel 1")
			library(tidyverse)
			library(ggplot2)
			library(ggbeeswarm)
			df_means <- exp_des %>% group_by(Label) %>% summarise(value=mean(value))

			pdf(paste("/public/ZhangJinLab/backdata/lahdata/TGene/",genename,"_CD80_boxplot.pdf",sep=""),width=8,height=10)
			p <- ggplot(data=exp_des, aes(x=Label,y=value)) + geom_quasirandom(aes(color=Label),method="frowney",width=0.3,size=1) + scale_color_manual(values=q4[c(1,2,4)]) + geom_boxplot(width=0.05,outlier.shape=NA,color="grey",fill="white",linetype=2) + stat_summary(fun.y=mean, geom="point", shape=20, size=10,  alpha=0.7, color="#8B814C")
			p <- p + xlab("") + ylab("Normalized expression level") + ggtitle("") + theme_classic()
			p <- p + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),panel.background=element_blank(), rect=element_rect(linetype=1), axis.line=element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11),legend.position="none") + theme(axis.text.x=element_text(vjust = 0.6, angle = 45))
			print(p)
			dev.off()
			iterer <- iterer + 1
		}
	}
}
write.table(res,file="/public/ZhangJinLab/backdata/lahdata/TGene/Tgenes_CD80_stat.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)


# Part 3
EnsembIDs <- read.table(file="/public/ZhangJinLab/backdata/lahdata/Sgenes",sep="\t",header=FALSE,stringsAsFactors=FALSE)
mapInfo <- read.table(file="/public/ZhangJinLab/backdata/lahdata/AA/features.tsv.gz",sep="\t",header=FALSE,stringsAsFactors=FALSE)
matchIndexes <- match(EnsembIDs[,2],mapInfo[,1])
geneVec <- mapInfo[matchIndexes[which(!is.na(matchIndexes))],2]
EnsembIDs <- EnsembIDs[which(!is.na(matchIndexes)),1]

res <- c()
iterer=1
for (genename in geneVec){
	if(length(which(rownames(normDataAA)==genename)) == 1){
		AA <- normDataAA[which(rownames(normDataAA)==genename),]
		AA <- AA[which(AA > 0)]
		cat(length(AA)/ncol(normDataAA),"\n")

		BB <- normDataBB[which(rownames(normDataBB)==genename),]
		BB <- BB[which(BB > 0)]
		cat(length(BB)/ncol(normDataBB),"\n")

		DD <- normDataDD[which(rownames(normDataDD)==genename),]
		DD <- DD[which(DD > 0)]
		cat(length(DD)/ncol(normDataDD),"\n")

		a <- wilcox.test(AA,BB)
		b <- wilcox.test(AA,DD)
		c <- wilcox.test(BB,DD)
		if(length(res) == 0){
			res <- c(genename,(length(AA)/ncol(normDataAA)),(length(BB)/ncol(normDataBB)),(length(DD)/ncol(normDataDD)),a$p.value,b$p.value,c$p.value)
		}else{
			res <- rbind(res,c(genename,(length(AA)/ncol(normDataAA)),(length(BB)/ncol(normDataBB)),(length(DD)/ncol(normDataDD)),a$p.value,b$p.value,c$p.value))
		}

		exp_des_label <- rep(c("AA","BB","DD"),c(length(AA),length(BB),length(DD)))
		exp_des_value <- as.numeric(c(AA,BB,DD))
		exp_des <- data.frame(cbind(exp_des_label,exp_des_value))
		colnames(exp_des) <- c("Label","Value")
		exp_des$Label=factor(exp_des_label,levels=c("AA","BB","DD"))
		exp_des$value <- as.numeric(exp_des_value)

		library(colorspace)
		q4 <- qualitative_hcl(4, palette = "Pastel 1")
		library(tidyverse)
		library(ggplot2)
		library(ggbeeswarm)
		df_means <- exp_des %>% group_by(Label) %>% summarise(value=mean(value))

		pdf(paste("/public/ZhangJinLab/backdata/lahdata/SGene/",genename,"_boxplot.pdf",sep=""),width=8,height=10)
		p <- ggplot(data=exp_des,aes(x=Label,y=value)) + geom_quasirandom(aes(color=Label),method="frowney",width=0.3,size=1) + scale_color_manual(values=q4[c(1,2,4)]) + geom_boxplot(width=0.05,outlier.shape=NA,color="grey",fill="white",linetype=2) + stat_summary(fun.y=mean, geom="point", shape=20, size=10,  alpha=0.7, color="#8B814C")
		p <- p + xlab("") + ylab("Normalized expression level") + ggtitle("") + theme_classic()
		p <- p + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),panel.background=element_blank(), rect=element_rect(linetype=1), axis.line=element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11),legend.position="none") + theme(axis.text.x=element_text(vjust = 0.6, angle = 45))
		print(p)
		dev.off()
		iterer <- iterer + 1
	}
}
write.table(res,file="/public/ZhangJinLab/backdata/lahdata/SGene/Sgenes_stat.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

res <- c()
iterer=1
for (genename in geneVec){
	if(length(which(rownames(normDataAA)==genename)) == 1){
		AA <- normDataAA[which(rownames(normDataAA)==genename),]
		CD80AA <- normDataAA["CD80",]
		AA <- AA[which(CD80AA > 0)]
		AA <- AA[which(AA > 0)]
		
		cat(length(AA)/ncol(normDataAA),"\n")

		BB <- normDataBB[which(rownames(normDataBB)==genename),]
		CD80BB <- normDataBB["CD80",]
		BB <- BB[which(CD80BB > 0)]
		BB <- BB[which(BB > 0)]
		cat(length(BB)/ncol(normDataBB),"\n")

		DD <- normDataDD[which(rownames(normDataDD)==genename),]
		CD80DD <- normDataDD["CD80",]
		DD <- DD[which(CD80DD > 0)]
		DD <- DD[which(DD > 0)]
		cat(length(DD)/ncol(normDataDD),"\n")
		
		if(length(AA) > 3 && length(BB) > 3 && length(DD) > 3){
			a <- wilcox.test(AA,BB)
			b <- wilcox.test(AA,DD)
			c <- wilcox.test(BB,DD)
			if(length(res) == 0){
				res <- c(genename,(length(AA)/ncol(normDataAA)),(length(BB)/ncol(normDataBB)),(length(DD)/ncol(normDataDD)),a$p.value,b$p.value,c$p.value)
			}else{
				res <- rbind(res,c(genename,(length(AA)/ncol(normDataAA)),(length(BB)/ncol(normDataBB)),(length(DD)/ncol(normDataDD)),a$p.value,b$p.value,c$p.value))
			}

			exp_des_label <- rep(c("AA","BB","DD"),c(length(AA),length(BB),length(DD)))
			exp_des_value <- as.numeric(c(AA,BB,DD))
			exp_des <- data.frame(cbind(exp_des_label,exp_des_value))
			colnames(exp_des) <- c("Label","Value")
			exp_des$Label=factor(exp_des_label,levels=c("AA","BB","DD"))
			exp_des$value <- as.numeric(exp_des_value)

			library(colorspace)
			q4 <- qualitative_hcl(4, palette = "Pastel 1")
			library(tidyverse)
			library(ggplot2)
			library(ggbeeswarm)
			df_means <- exp_des %>% group_by(Label) %>% summarise(value=mean(value))

			pdf(paste("/public/ZhangJinLab/backdata/lahdata/SGene/",genename,"_CD80_boxplot.pdf",sep=""),width=8,height=10)
			p <- ggplot(data=exp_des, aes(x=Label,y=value)) + geom_quasirandom(aes(color=Label),method="frowney",width=0.3,size=1) + scale_color_manual(values=q4[c(1,2,4)]) + geom_boxplot(width=0.05,outlier.shape=NA,color="grey",fill="white",linetype=2) + stat_summary(fun.y=mean, geom="point", shape=20, size=10,  alpha=0.7, color="#8B814C")
			p <- p + xlab("") + ylab("Normalized expression level") + ggtitle("") + theme_classic()
			p <- p + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),panel.background=element_blank(), rect=element_rect(linetype=1), axis.line=element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11),legend.position="none") + theme(axis.text.x=element_text(vjust = 0.6, angle = 45))
			print(p)
			dev.off()
			iterer <- iterer + 1
		}
	}
}
write.table(res,file="/public/ZhangJinLab/backdata/lahdata/SGene/Sgenes_CD80_stat.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
