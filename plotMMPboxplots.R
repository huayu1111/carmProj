library(Seurat)
library(ggplot2)
library(beeswarm)
load("/media/data/AA/outs/filtered_feature_bc_matrix/10X_LAH_InHouse.RData")
normDataAA <-  GetAssayData(object = mydata_total[["RNA"]], slot = "data")
load("/media/data/BB/outs/filtered_feature_bc_matrix/10X_LAH_InHouse.RData")
normDataBB <-  GetAssayData(object = mydata_total[["RNA"]], slot = "data")
load("/media/data/DD/outs/filtered_feature_bc_matrix/10X_LAH_InHouse.RData")
normDataDD <-  GetAssayData(object = mydata_total[["RNA"]], slot = "data")
geneVec <- c("MMP1","MMP10","MMP11","MMP12","MMP13","MMP14","MMP15","MMP16","MMP17","MMP19","MMP2","MMP20","MMP21","MMP23A","MMP23B","MMP24","MMP25","MMP26","MMP27","MMP28","MMP3","MMP7","MMP8","MMP9","TIMP1","TIMP2","TIMP3","TIMP4")

res <- c()
for (genename in geneVec){
	if(length(which(rownames(normDataAA)==genename)) == 1){
		CXCL9AA <- normDataAA[which(rownames(normDataAA)==genename),]
		CXCL9AA <- CXCL9AA[which(CXCL9AA > 0)]
		cat(length(CXCL9AA)/ncol(normDataAA),"\n")

		CXCL9BB <- normDataBB[which(rownames(normDataBB)==genename),]
		CXCL9BB <- CXCL9BB[which(CXCL9BB > 0)]
		cat(length(CXCL9BB)/ncol(normDataBB),"\n")

		CXCL9DD <- normDataDD[which(rownames(normDataDD)==genename),]
		CXCL9DD <- CXCL9DD[which(CXCL9DD > 0)]
		cat(length(CXCL9DD)/ncol(normDataDD),"\n")

		a <- wilcox.test(CXCL9AA,CXCL9BB)
		b <- wilcox.test(CXCL9AA,CXCL9DD)
		c <- wilcox.test(CXCL9BB,CXCL9DD)
		if(length(res) == 0){
			res <- c(genename,(length(CXCL9AA)/ncol(normDataAA)),(length(CXCL9BB)/ncol(normDataBB)),(length(CXCL9DD)/ncol(normDataDD)),a$p.value,b$p.value,c$p.value)
		}else{
			res <- rbind(res,c(genename,(length(CXCL9AA)/ncol(normDataAA)),(length(CXCL9BB)/ncol(normDataBB)),(length(CXCL9DD)/ncol(normDataDD)),a$p.value,b$p.value,c$p.value))
		}

		exp_des_label <- rep(c("AA","BB","DD"),c(length(CXCL9AA),length(CXCL9BB),length(CXCL9DD)))
		exp_des_value <- as.numeric(c(CXCL9AA,CXCL9BB,CXCL9DD))
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

		pdf(paste("/media/data/AA/",genename,"_boxplot.pdf",sep=""),width=8,height=10)
		p <- ggplot(data=exp_des, aes(x=Label,y=value)) + geom_quasirandom(aes(color=Label),method="frowney",width=0.3,size=1) + scale_color_manual(values=q4[c(1,2,4)]) + geom_boxplot(width=0.05,outlier.shape=NA,color="grey",fill="white",linetype=2) + stat_summary(fun.y=mean, geom="point", shape=20, size=10,  alpha=0.7, color="#8B814C")
		p <- p + xlab("") + ylab("Normalized expression level") + ggtitle("") + theme_classic()
		p <- p + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),panel.background=element_blank(), rect=element_rect(linetype=1), axis.line=element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11),legend.position="none") + theme(axis.text.x=element_text(vjust = 0.6, angle = 45))
		print(p)
		dev.off()
	}
}
write.table(res,file="/media/data/AA/MMP_stat.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

res <- c()
for (genename in geneVec){
	if(length(which(rownames(normDataAA)==genename)) == 1){
		CXCL9AA <- normDataAA[which(rownames(normDataAA)==genename),]
		CD80AA <- normDataAA["CD80",]
		CXCL9AA <- CXCL9AA[which(CD80AA > 0)]
		CXCL9AA <- CXCL9AA[which(CXCL9AA > 0)]
		
		cat(length(CXCL9AA)/ncol(normDataAA),"\n")

		CXCL9BB <- normDataBB[which(rownames(normDataBB)==genename),]
		CD80BB <- normDataBB["CD80",]
		CXCL9BB <- CXCL9BB[which(CD80BB > 0)]
		CXCL9BB <- CXCL9BB[which(CXCL9BB > 0)]
		cat(length(CXCL9BB)/ncol(normDataBB),"\n")

		CXCL9DD <- normDataDD[which(rownames(normDataDD)==genename),]
		CD80DD <- normDataDD["CD80",]
		CXCL9DD <- CXCL9DD[which(CD80DD > 0)]
		CXCL9DD <- CXCL9DD[which(CXCL9DD > 0)]
		cat(length(CXCL9DD)/ncol(normDataDD),"\n")
		
		if(length(CXCL9AA) > 3 && length(CXCL9BB) > 3 && length(CXCL9DD) > 3){
			a <- wilcox.test(CXCL9AA,CXCL9BB)
			b <- wilcox.test(CXCL9AA,CXCL9DD)
			c <- wilcox.test(CXCL9BB,CXCL9DD)
			if(length(res) == 0){
				res <- c(genename,(length(CXCL9AA)/ncol(normDataAA)),(length(CXCL9BB)/ncol(normDataBB)),(length(CXCL9DD)/ncol(normDataDD)),a$p.value,b$p.value,c$p.value)
			}else{
				res <- rbind(res,c(genename,(length(CXCL9AA)/ncol(normDataAA)),(length(CXCL9BB)/ncol(normDataBB)),(length(CXCL9DD)/ncol(normDataDD)),a$p.value,b$p.value,c$p.value))
			}

			exp_des_label <- rep(c("AA","BB","DD"),c(length(CXCL9AA),length(CXCL9BB),length(CXCL9DD)))
			exp_des_value <- as.numeric(c(CXCL9AA,CXCL9BB,CXCL9DD))
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

			pdf(paste("/media/data/AA/",genename,"_CD80_boxplot.pdf",sep=""),width=8,height=10)
			p <- ggplot(data=exp_des, aes(x=Label,y=value)) + geom_quasirandom(aes(color=Label),method="frowney",width=0.3,size=1) + scale_color_manual(values=q4[c(1,2,4)]) + geom_boxplot(width=0.05,outlier.shape=NA,color="grey",fill="white",linetype=2) + stat_summary(fun.y=mean, geom="point", shape=20, size=10,  alpha=0.7, color="#8B814C")
			p <- p + xlab("") + ylab("Normalized expression level") + ggtitle("") + theme_classic()
			p <- p + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),panel.background=element_blank(), rect=element_rect(linetype=1), axis.line=element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11),legend.position="none") + theme(axis.text.x=element_text(vjust = 0.6, angle = 45))
			print(p)
			dev.off()
		}
	}
}
write.table(res,file="/media/data/AA/MMP_CD80_stat.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
