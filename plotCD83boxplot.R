library(Seurat)
library(ggplot2)
library(beeswarm)

load("/media/data/AA/outs/filtered_feature_bc_matrix/10X_LAH_InHouse.RData")
normDataAA <-  GetAssayData(object = mydata_total[["RNA"]], slot = "data")
CD83AA <- normDataAA["CD83",]
CD83AA <- CD83AA[which(CD83AA > 0)]
cat(length(CD83AA)/ncol(normDataAA))

load("/media/data/BB/outs/filtered_feature_bc_matrix/10X_LAH_InHouse.RData")
normDataBB <-  GetAssayData(object = mydata_total[["RNA"]], slot = "data")
CD83BB <- normDataBB["CD83",]
CD83BB <- CD83BB[which(CD83BB > 0)]
cat(length(CD83BB)/ncol(normDataBB))

load("/media/data/DD/outs/filtered_feature_bc_matrix/10X_LAH_InHouse.RData")
normDataDD <-  GetAssayData(object = mydata_total[["RNA"]], slot = "data")
CD83DD <- normDataDD["CD83",]
CD83DD <- CD83DD[which(CD83DD > 0)]
cat(length(CD83DD)/ncol(normDataDD))

wilcox.test(CD83AA,CD83BB)
wilcox.test(CD83AA,CD83DD)
wilcox.test(CD83BB,CD83DD)

exp_des_label <- rep(c("AA","BB","DD"),c(length(CD83AA),length(CD83BB),length(CD83DD)))
exp_des_value <- as.numeric(c(CD83AA,CD83BB,CD83DD))
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

pdf("/media/data/AA/CD83_boxplot.pdf",width=8,height=10)
p <- ggplot(data=exp_des, aes(x=Label,y=value)) + geom_quasirandom(aes(color=Label),method="frowney",width=0.3,size=1) + scale_color_manual(values=q4[c(1,2,4)]) + geom_boxplot(width=0.05,outlier.shape=NA,color="grey",fill="white",linetype=2) + stat_summary(fun.y=mean, geom="point", shape=20, size=10,  alpha=0.7, color="#8B814C")
p <- p + xlab("") + ylab("Normalized expression level") + ggtitle("") + theme_classic()
p <- p + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),panel.background=element_blank(), rect=element_rect(linetype=1), axis.line=element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11),legend.position="none") + theme(axis.text.x=element_text(vjust = 0.6, angle = 45))
print(p)
dev.off()


library(Seurat)
library(ggplot2)
library(beeswarm)

load("/media/data/AA/outs/filtered_feature_bc_matrix/10X_LAH_InHouse.RData")
normDataAA <-  GetAssayData(object = mydata_total[["RNA"]], slot = "data")
CD83AA <- normDataAA["CD83",]
CD80AA <- normDataAA["CD80",]
CD83AA <- CD83AA[which(CD80AA > 0)]
CD83AA <- CD83AA[which(CD83AA > 0)]

load("/media/data/BB/outs/filtered_feature_bc_matrix/10X_LAH_InHouse.RData")
normDataBB <-  GetAssayData(object = mydata_total[["RNA"]], slot = "data")
CD83BB <- normDataBB["CD83",]
CD80BB <- normDataBB["CD80",]
CD83BB <- CD83BB[which(CD80BB > 0)]
CD83BB <- CD83BB[which(CD83BB > 0)]

load("/media/data/DD/outs/filtered_feature_bc_matrix/10X_LAH_InHouse.RData")
normDataDD <-  GetAssayData(object = mydata_total[["RNA"]], slot = "data")
CD83DD <- normDataDD["CD83",]
CD80DD <- normDataDD["CD80",]
CD83DD <- CD83DD[which(CD80DD > 0)]
CD83DD <- CD83DD[which(CD83DD > 0)]

wilcox.test(CD83AA,CD83BB)
wilcox.test(CD83AA,CD83DD)
wilcox.test(CD83BB,CD83DD)

exp_des_label <- rep(c("AA","BB","DD"),c(length(CD83AA),length(CD83BB),length(CD83DD)))
exp_des_value <- as.numeric(c(CD83AA,CD83BB,CD83DD))
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

pdf("/media/data/AA/CD83_CD80Pos_boxplot.pdf",width=8,height=10)
p <- ggplot(data=exp_des, aes(x=Label,y=value)) + geom_quasirandom(aes(color=Label),method="frowney",width=0.3,size=1) + scale_color_manual(values=q4[c(1,2,4)]) + geom_boxplot(width=0.05,outlier.shape=NA,color="grey",fill="white",linetype=2) + stat_summary(fun.y=mean, geom="point", shape=20, size=10,  alpha=0.7, color="#8B814C")
p <- p + xlab("") + ylab("Normalized expression level") + ggtitle("") + theme_classic()
p <- p + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),panel.background=element_blank(), rect=element_rect(linetype=1), axis.line=element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11),legend.position="none") + theme(axis.text.x=element_text(vjust = 0.6, angle = 45))
print(p)
dev.off()

