library(Seurat)
library(ggplot2)
library(beeswarm)

load("/media/data/AA/outs/filtered_feature_bc_matrix/10X_LAH_InHouse.RData")
normDataAA <-  GetAssayData(object = mydata_total[["RNA"]], slot = "data")
CD80AA <- normDataAA["CD80",]
CD80AA <- CD80AA[which(CD80AA > 0)]
cat(length(CD80AA)/ncol(normDataAA))

load("/media/data/BB/outs/filtered_feature_bc_matrix/10X_LAH_InHouse.RData")
normDataBB <-  GetAssayData(object = mydata_total[["RNA"]], slot = "data")
CD80BB <- normDataBB["CD80",]
CD80BB <- CD80BB[which(CD80BB > 0)]
cat(length(CD80BB)/ncol(normDataBB))

load("/media/data/DD/outs/filtered_feature_bc_matrix/10X_LAH_InHouse.RData")
normDataDD <-  GetAssayData(object = mydata_total[["RNA"]], slot = "data")
CD80DD <- normDataDD["CD80",]
CD80DD <- CD80DD[which(CD80DD > 0)]
cat(length(CD80DD)/ncol(normDataDD))

wilcox.test(CD80AA,CD80BB)
wilcox.test(CD80AA,CD80DD)
wilcox.test(CD80BB,CD80DD)

exp_des_label <- rep(c("AA","BB","DD"),c(length(CD80AA),length(CD80BB),length(CD80DD)))
exp_des_value <- as.numeric(c(CD80AA,CD80BB,CD80DD))
exp_des <- data.frame(cbind(exp_des_label,exp_des_value))
colnames(exp_des) <- c("Label","Value")
exp_des$Label=factor(exp_des_label,levels=c("AA","BB","DD"))
exp_des$value <- as.numeric(exp_des_value)

pdf("/media/data/AA/CD80_boxplot.pdf",width=8,height=10)
boxplot(exp_des$value ~ as.factor(exp_des$Label),outline = FALSE,xlab="",ylab="Normalized expression level")
beeswarm(exp_des$value ~ as.factor(exp_des$Label),xlab="",ylab="Normalized expression level",col=c("darkorchid3","goldenrod1","dodgerblue4"),pch=16,cex=0.3,add=TRUE)
dev.off()
