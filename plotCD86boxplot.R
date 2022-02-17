library(Seurat)
library(ggplot2)
library(beeswarm)

load("/media/data/AA/outs/filtered_feature_bc_matrix/10X_LAH_InHouse.RData")
normDataAA <-  GetAssayData(object = mydata_total[["RNA"]], slot = "data")
CD86AA <- normDataAA["CD86",]
CD86AA <- CD86AA[which(CD86AA > 0)]
cat(length(CD86AA)/ncol(normDataAA))

load("/media/data/BB/outs/filtered_feature_bc_matrix/10X_LAH_InHouse.RData")
normDataBB <-  GetAssayData(object = mydata_total[["RNA"]], slot = "data")
CD86BB <- normDataBB["CD86",]
CD86BB <- CD86BB[which(CD86BB > 0)]
cat(length(CD86BB)/ncol(normDataBB))

load("/media/data/DD/outs/filtered_feature_bc_matrix/10X_LAH_InHouse.RData")
normDataDD <-  GetAssayData(object = mydata_total[["RNA"]], slot = "data")
CD86DD <- normDataDD["CD86",]
CD86DD <- CD86DD[which(CD86DD > 0)]
cat(length(CD86DD)/ncol(normDataDD))

wilcox.test(CD86AA,CD86BB)
wilcox.test(CD86AA,CD86DD)
wilcox.test(CD86BB,CD86DD)

exp_des_label <- rep(c("AA","BB","DD"),c(length(CD86AA),length(CD86BB),length(CD86DD)))
exp_des_value <- as.numeric(c(CD86AA,CD86BB,CD86DD))
exp_des <- data.frame(cbind(exp_des_label,exp_des_value))
colnames(exp_des) <- c("Label","Value")
exp_des$Label=factor(exp_des_label,levels=c("AA","BB","DD"))
exp_des$value <- as.numeric(exp_des_value)

pdf("/media/data/AA/CD86_boxplot.pdf",width=8,height=10)
boxplot(exp_des$value ~ as.factor(exp_des$Label),outline = FALSE,xlab="",ylab="Normalized expression level")
beeswarm(exp_des$value ~ as.factor(exp_des$Label),xlab="",ylab="Normalized expression level",col=c("darkorchid3","goldenrod1","dodgerblue4"),pch=16,cex=0.3,add=TRUE)
dev.off()


