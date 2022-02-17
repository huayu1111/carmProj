library(Matrix)
library(Seurat)
library(dplyr)
library(reshape2)
library(ggplot2)

norm_cross_col <- function(vec_value){
	vec_value_return <- (vec_value-mean(vec_value))/sd(vec_value)
}

MarkerInfo_HMDM <- read.delim("/media/data/lahannodata/HMDM_M1M2_markers.txt", header = TRUE, stringsAsFactors = FALSE)
MarkerInfo_IPSDM <- read.delim("/media/data/lahannodata/IPSDM_M1M2_markers.txt", header = TRUE, stringsAsFactors = FALSE)

MarkerInfo <- cbind(MarkerInfo_HMDM,MarkerInfo_IPSDM[,seq(3,12)])
M1M2Expr <- t(apply(MarkerInfo[,seq(3,14)],1,norm_cross_col))
M1M2order <- order(M1M2Expr[,1],decreasing=TRUE)
G1 <- MarkerInfo[M1M2order,2]

mapInfo <- read.table(file = "/media/data/lahannodata/gene_TPM.txt", sep = "\t", header = T, row.names = 1)
MarkerInfo <- read.delim("/media/data/lahannodata/diffanaly_markers.txt", header = TRUE, stringsAsFactors = FALSE)
matchIndexes <- match(MarkerInfo[,2],mapInfo[,2])
MarkerInfo[,seq(3,ncol(MarkerInfo))] <- mapInfo[matchIndexes,c("SRR1182376","SRR1182388","SRR2939146","SRR1182390","SRR2910670","SRR2910671","SRR1182378","SRR1182392","SRR1182394","SRR2939152","SRR2910672","SRR2910673")]
M1M2Expr <- t(apply(MarkerInfo[,seq(3,14)],1,norm_cross_col))
M1M2order <- order(M1M2Expr[,1],decreasing=TRUE)
G2 <- MarkerInfo[M1M2order,2]


mydata_total_dir <- "/media/data/DD/outs/filtered_feature_bc_matrix/"
mydata.barcode.path <- paste0(mydata_total_dir, "barcodes.tsv.gz")
mydata.features.path <- paste0(mydata_total_dir, "features.tsv.gz")
mydata.matrix.path <- paste0(mydata_total_dir, "matrix.mtx.gz")
mydata_total <- readMM(file = mydata.matrix.path)
mydata.feature.names = read.delim(mydata.features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
mydata.barcode.names = read.delim(mydata.barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mydata_total) = mydata.barcode.names$V1
rownames(mydata_total) = mydata.feature.names$V2

mydata_total <- CreateSeuratObject(counts = mydata_total, project = "CART_10X_LX5", min.cells = 3,min.features = 200)
mydata_total[["percent.mt"]] <- PercentageFeatureSet(mydata_total, pattern = "^MT-")
mydata_total <- subset(mydata_total, subset = nFeature_RNA > 200 & nFeature_RNA < 30000 & percent.mt < 10)
mydata_total <- NormalizeData(mydata_total, normalization.method = "LogNormalize", scale.factor = 10000)
# mydata_total <- FindVariableFeatures(mydata_total, selection.method = "vst", nfeatures = 5000)
var.genes <- c(G1,G2)
mydata_total <- ScaleData(mydata_total, vars.to.regress = c("nFeature_RNA", "percent.mt"), display.progress = T, do.par=TRUE, num.cores=40)
mydata_total <- RunPCA(mydata_total, features = var.genes)
mydata_total <- RunTSNE(mydata_total, dims = 1:10)
mydata_total <- RunUMAP(mydata_total, dims = 1:10)
mydata_total <- FindNeighbors(mydata_total,dims = 1:10)
mydata_total <- FindClusters(mydata_total, resolution = 0.25)
save(mydata_total, file="/media/data/DD/outs/filtered_feature_bc_matrix/10X_LAH_InHouse.RData")

pdf("/media/data/DD/outs/filtered_feature_bc_matrix/Scatter_UMAP_dim1VSdim2_InHouse.pdf",width=10,height=10)
DimPlot(mydata_total, reduction = 'umap', pt.size = 0.5)
dev.off()

pdf("/media/data/DD/outs/filtered_feature_bc_matrix/Scatter_tSNE_dim1VSdim2_InHouse.pdf",width=10,height=10)
DimPlot(mydata_total, reduction = 'tsne', pt.size = 0.5)
dev.off()

for (genename in var.genes){
	FeaturePlot(object = mydata_total,features = genename,cols = c("grey", "red"),reduction = "umap")
	ggsave(paste("/media/data/DD/outs/filtered_feature_bc_matrix/",genename,".pdf",sep=""),width=5,height=5)
}
