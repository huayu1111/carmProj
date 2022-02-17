library(Seurat)
library(ggplot2)
library(beeswarm)

load("/media/data/AA/outs/filtered_feature_bc_matrix/10X_LAH_InHouse.RData")
normDataAA <-  GetAssayData(object = mydata_total[["RNA"]], slot = "data")
CD206AA <- normDataAA["MRC1",]
CD206AA <- CD206AA[which(CD206AA > 0)]
cat(length(CD206AA)/ncol(normDataAA))

load("/media/data/BB/outs/filtered_feature_bc_matrix/10X_LAH_InHouse.RData")
normDataBB <-  GetAssayData(object = mydata_total[["RNA"]], slot = "data")
CD206BB <- normDataBB["MRC1",]
CD206BB <- CD206BB[which(CD206BB > 0)]
cat(length(CD206BB)/ncol(normDataBB))
load("/media/data/DD/outs/filtered_feature_bc_matrix/10X_LAH_InHouse.RData")
normDataDD <-  GetAssayData(object = mydata_total[["RNA"]], slot = "data")
CD206DD <- normDataDD["MRC1",]
CD206DD <- CD206DD[which(CD206DD > 0)]
cat(length(CD206DD)/ncol(normDataDD))

wilcox.test(CD206AA,CD206BB)
wilcox.test(CD206AA,CD206DD)
wilcox.test(CD206BB,CD206DD)

exp_des_label <- rep(c("AA","BB","DD"),c(length(CD206AA),length(CD206BB),length(CD206DD)))
exp_des_value <- as.numeric(c(CD206AA,CD206BB,CD206DD))
exp_des <- data.frame(cbind(exp_des_label,exp_des_value))
colnames(exp_des) <- c("Label","Value")
exp_des$Label=factor(exp_des_label,levels=c("AA","BB","DD"))
exp_des$value <- as.numeric(exp_des_value)

pdf("/media/data/AA/CD206_boxplot.pdf",width=8,height=10)
boxplot(exp_des$value ~ as.factor(exp_des$Label),outline = FALSE,xlab="",ylab="Normalized expression level")
beeswarm(exp_des$value ~ as.factor(exp_des$Label),xlab="",ylab="Normalized expression level",col=c("darkorchid3","goldenrod1","dodgerblue4"),pch=16,cex=0.3,add=TRUE)
dev.off()
