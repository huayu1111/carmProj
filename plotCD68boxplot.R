library(Seurat)
library(ggplot2)
library(beeswarm)

load("/media/data/AA/outs/filtered_feature_bc_matrix/10X_LAH_InHouse.RData")
normDataAA <-  GetAssayData(object = mydata_total[["RNA"]], slot = "data")
CD68AA <- normDataAA["CD68",]
CD68AA <- CD68AA[which(CD68AA > 0)]

load("/media/data/BB/outs/filtered_feature_bc_matrix/10X_LAH_InHouse.RData")
normDataBB <-  GetAssayData(object = mydata_total[["RNA"]], slot = "data")
CD68BB <- normDataBB["CD68",]
CD68BB <- CD68BB[which(CD68BB > 0)]

load("/media/data/DD/outs/filtered_feature_bc_matrix/10X_LAH_InHouse.RData")
normDataDD <-  GetAssayData(object = mydata_total[["RNA"]], slot = "data")
CD68DD <- normDataDD["CD68",]
CD68DD <- CD68DD[which(CD68DD > 0)]

wilcox.test(CD68AA,CD68BB)
wilcox.test(CD68AA,CD68DD)
wilcox.test(CD68BB,CD68DD)

exp_des_label <- rep(c("AA","BB","DD"),c(length(CD68AA),length(CD68BB),length(CD68DD)))
exp_des_value <- as.numeric(c(CD68AA,CD68BB,CD68DD))
exp_des <- data.frame(cbind(exp_des_label,exp_des_value))
colnames(exp_des) <- c("Label","Value")
exp_des$Label=factor(exp_des_label,levels=c("AA","BB","DD"))
exp_des$value <- as.numeric(exp_des_value)

pdf("/media/data/AA/CD68_boxplot.pdf",width=8,height=10)
boxplot(exp_des$value ~ as.factor(exp_des$Label),outline = FALSE,xlab="",ylab="Normalized expression level")
beeswarm(exp_des$value ~ as.factor(exp_des$Label),xlab="",ylab="Normalized expression level",col=c("darkorchid3","goldenrod1","dodgerblue4"),pch=16,cex=0.3,add=TRUE)
dev.off()