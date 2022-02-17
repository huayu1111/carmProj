library(reshape2)
A <- read.table(file="/media/data/AA/M1Perc.txt",header=TRUE,row.names=1)
A <- as.matrix(A)
A <- A[,c(1,3,4,5,6)]
df <- data.frame(C=rep(row.names(A),ncol(A)),G=rep(colnames(A),each=nrow(A)),P=as.vector(A))

df$G <- factor(df$G,levels=c("CD80","CD83","CCL8","CXCL9","CXCL11")) 
df$C <- factor(df$C,levels=rev(c("AA","BB","DD")))
df$P <- as.numeric(A)

library(ggplot2)
pdf("/media/data/AA/M1marker_perc.pdf",width=10,height=5)
p <- ggplot(data=df, aes(x=C, y=P, fill=G)) + geom_bar(stat="identity",position=position_dodge())
p <- p + theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) 
p <- p + coord_flip()
print(p)
dev.off()


library(reshape2)
A <- read.table(file="/media/data/AA/M2Perc.txt",header=TRUE,row.names=1)
A <- as.matrix(A)
df <- data.frame(C=rep(row.names(A),ncol(A)),G=rep(colnames(A),each=nrow(A)),P=as.vector(A))

df$G <- factor(df$G,levels=c("CD206","CCL13","CSF1R"))
df$C <- factor(df$C,levels=rev(c("AA","BB","DD")))
df$P <- as.numeric(A)

library(ggplot2)
pdf("/media/data/AA/M2marker_perc.pdf",width=10,height=5)
p <- ggplot(data=df, aes(x=C, y=P, fill=G)) + geom_bar(stat="identity",position=position_dodge())
p <- p + theme(panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) 
p <- p + coord_flip()
print(p)
dev.off()