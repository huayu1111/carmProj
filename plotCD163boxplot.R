library(Seurat)
library(ggplot2)
library(beeswarm)

load("/media/data/AA/outs/filtered_feature_bc_matrix/10X_LAH_InHouse.RData")
normDataAA <-  GetAssayData(object = mydata_total[["RNA"]], slot = "data")
CD163AA <- normDataAA["CD163",]
# CD163AA <- CD163AA[which(CD163AA > 0)]
# cat(length(CD163AA)/ncol(normDataAA))

load("/media/data/BB/outs/filtered_feature_bc_matrix/10X_LAH_InHouse.RData")
normDataBB <-  GetAssayData(object = mydata_total[["RNA"]], slot = "data")
CD163BB <- normDataBB["CD163",]
# CD163BB <- CD163BB[which(CD163BB > 0)]
# cat(length(CD163BB)/ncol(normDataBB))

load("/media/data/DD/outs/filtered_feature_bc_matrix/10X_LAH_InHouse.RData")
normDataDD <-  GetAssayData(object = mydata_total[["RNA"]], slot = "data")
CD163DD <- normDataDD["CD163",]
# CD163DD <- CD163DD[which(CD163DD > 0)]
# cat(length(CD163DD)/ncol(normDataDD))

wilcox.test(CD163AA,CD163BB)
wilcox.test(CD163AA,CD163DD)
wilcox.test(CD163BB,CD163DD)

exp_des_label <- rep(c("AA","BB","DD"),c(length(CD163AA),length(CD163BB),length(CD163DD)))
exp_des_value <- as.numeric(c(CD163AA,CD163BB,CD163DD))
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

pdf("/media/data/AA/CD163_boxplot.pdf",width=8,height=10)
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
CD163AA <- normDataAA["CD163",]
MRC1AA <- normDataAA["MRC1",]
CD163AA <- CD163AA[which(MRC1AA > 0)]
# CD163AA <- CD163AA[which(CD163AA > 0)]

load("/media/data/BB/outs/filtered_feature_bc_matrix/10X_LAH_InHouse.RData")
normDataBB <-  GetAssayData(object = mydata_total[["RNA"]], slot = "data")
CD163BB <- normDataBB["CD163",]
MRC1BB <- normDataBB["MRC1",]
CD163BB <- CD163BB[which(MRC1BB > 0)]
# CD163BB <- CD163BB[which(CD163BB > 0)]

load("/media/data/DD/outs/filtered_feature_bc_matrix/10X_LAH_InHouse.RData")
normDataDD <-  GetAssayData(object = mydata_total[["RNA"]], slot = "data")
CD163DD <- normDataDD["CD163",]
MRC1DD <- normDataDD["MRC1",]
CD163DD <- CD163DD[which(MRC1DD > 0)]
# CD163DD <- CD163DD[which(CD163DD > 0)]

wilcox.test(CD163AA,CD163BB)
wilcox.test(CD163AA,CD163DD)
wilcox.test(CD163BB,CD163DD)

exp_des_label <- rep(c("AA","BB","DD"),c(length(CD163AA),length(CD163BB),length(CD163DD)))
exp_des_value <- as.numeric(c(CD163AA,CD163BB,CD163DD))
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

pdf("/media/data/AA/CD163_MRC1Pos_boxplot.pdf",width=8,height=10)
p <- ggplot(data=exp_des, aes(x=Label,y=value)) + geom_quasirandom(aes(color=Label),method="frowney",width=0.3,size=1) + scale_color_manual(values=q4[c(1,2,4)]) + geom_boxplot(width=0.05,outlier.shape=NA,color="grey",fill="white",linetype=2) + stat_summary(fun.y=mean, geom="point", shape=20, size=10,  alpha=0.7, color="#8B814C")
p <- p + xlab("") + ylab("Normalized expression level") + ggtitle("") + theme_classic()
p <- p + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(),panel.background=element_blank(), rect=element_rect(linetype=1), axis.line=element_line(colour = "black"),axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=11),axis.text.y=element_text(size=11),legend.position="none") + theme(axis.text.x=element_text(vjust = 0.6, angle = 45))
print(p)
dev.off()

