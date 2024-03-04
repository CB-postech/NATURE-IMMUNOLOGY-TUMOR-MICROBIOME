library(SingleCellExperiment)
library(DropletUtils)
library(biomaRt)
library(scater)
library(scran)
library(biomaRt)
library(RColorBrewer)
library(Seurat)

setwd("D:/Postech_TIL/data/")

dir_name_1 <- "D:/Postech_TIL/data/raw/20190906_C10/outs/raw_feature_bc_matrix/"
dir_name_2 <- "D:/Postech_TIL/data/raw/20190906_C11/outs/raw_feature_bc_matrix/"
dir_name_3 <- "D:/Postech_TIL/data/raw/20190906_C12/outs/raw_feature_bc_matrix/"


C10_sce <- read10xCounts(dir_name_1)
C11_sce <- read10xCounts(dir_name_2)
C12_sce <- read10xCounts(dir_name_3)




# QC for C10
br.out <- barcodeRanks(counts(C10_sce))

plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")

set.seed(100)
e.out <- emptyDrops(counts(C10_sce))
e.out
is.cell <- e.out$FDR <= 0.05
sum(is.cell, na.rm=TRUE)

abline(h=br.out[br.out$rank == sum(is.cell, na.rm=TRUE)+1,]$total, col="red", lty=2)
abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen", "red"), 
       legend=c("knee", "inflection", "FDR_0.05"))


C10_sce <- C10_sce[,which(e.out$FDR <= 0.05)]
C10_sce@assays$data$logcounts <- log2(C10_sce@assays$data$counts + 1)


ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl", host="www.ensembl.org")
ensemblGenes <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name',  'chromosome_name', 'gene_biotype'), mart=ensembl)
rownames(ensemblGenes) <- ensemblGenes[,1]

mtGenes <- ensemblGenes[ensemblGenes[,3]=="MT",]
is.mito <- rownames(C10_sce) %in% mtGenes[,1]
length(is.mito[is.mito == TRUE])
C10_sce <- calculateQCMetrics(C10_sce, feature_controls=list(Mito=is.mito))


h <- hist(C10_sce$log10_total_counts, 
          breaks=100, col="grey80",
          xlab="Log-total UMI count")
cuts <- cut(h$breaks, c(0,3.0,Inf))
plot(h, col=c("red","white")[cuts])

h <- hist(C10_sce$pct_counts_Mito, 
          breaks=100, col="grey80",
          xlab="Proportion of reads in mitochondrial genes")
cuts <- cut(h$breaks, c(10,100,Inf))
plot(h, col=c("red","white")[cuts])



C10_sce <- runPCA(C10_sce, use_coldata=TRUE)
plotReducedDim(C10_sce, use_dimred = "PCA_coldata")

ggplot(as.data.frame(C10_sce@reducedDims$PCA_coldata), 
       aes(x=PC1, y=PC2, color = C10_sce$pct_counts_Mito)) +
   scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) + 
   geom_point(size = 1) + 
   theme_bw() + 
   theme(text = element_text(size = 20),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         axis.line = element_line(size = 1),
         axis.ticks = element_line(size = 1),
         legend.title = element_blank(),
         legend.key = element_blank())

ggplot(as.data.frame(C10_sce@reducedDims$PCA_coldata), 
       aes(x=PC1, y=PC2, color = C10_sce$log10_total_counts)) +
   scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) + 
   geom_point(size = 0.1) + 
   theme_bw() + 
   theme(text = element_text(size = 20),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         axis.line = element_line(size = 1),
         axis.ticks = element_line(size = 1),
         legend.title = element_blank(),
         legend.key = element_blank())


C10_sce$use <- (
   (C10_sce$pct_counts_Mito <= 10) &
      (C10_sce$log10_total_counts >= 3)
)
table(C10_sce$use)

ggplot(as.data.frame(C10_sce@reducedDims$PCA_coldata), 
       aes(x=PC1, y=PC2, color = as.factor(C10_sce$use))) +
   geom_point(size = 1) + 
   theme_bw() + 
   theme(text = element_text(size = 20),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         axis.line = element_line(size = 1),
         axis.ticks = element_line(size = 1),
         legend.title = element_blank(),
         legend.key = element_blank())

sce_C10_filtered = C10_sce[,C10_sce$use]
save(sce_C10_filtered, file = "C10_sce_filtered.RData")




# QC for C11
br.out <- barcodeRanks(counts(C11_sce))

plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")

set.seed(100)
e.out <- emptyDrops(counts(C11_sce))
e.out
is.cell <- e.out$FDR <= 0.05
sum(is.cell, na.rm=TRUE)

abline(h=br.out[br.out$rank == sum(is.cell, na.rm=TRUE),]$total-0.5, col="red", lty=2)
abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen", "red"), 
       legend=c("knee", "inflection", "FDR_0.05"))


C11_sce <- C11_sce[,which(e.out$FDR <= 0.05)]
C11_sce@assays$data$logcounts <- log2(C11_sce@assays$data$counts + 1)


is.mito = rownames(C11_sce) %in% mtGenes[,1]
length(is.mito[is.mito == TRUE])
C11_sce = calculateQCMetrics(C11_sce, feature_controls=list(Mito=is.mito))


h <- hist(C11_sce$log10_total_counts, 
          breaks=100, col="grey80",
          xlab="Log-total UMI count")
cuts <- cut(h$breaks, c(0,3.0,Inf))
plot(h, col=c("red","white")[cuts])

h <- hist(C11_sce$pct_counts_Mito, 
          breaks=100, col="grey80",
          xlab="Proportion of reads in mitochondrial genes")
cuts <- cut(h$breaks, c(10,100,Inf))
plot(h, col=c("red","white")[cuts])


C11_sce <- runPCA(C11_sce, use_coldata=TRUE)
plotReducedDim(C11_sce, use_dimred = "PCA_coldata")

ggplot(as.data.frame(C11_sce@reducedDims$PCA_coldata), 
       aes(x=PC1, y=PC2, color = C11_sce$pct_counts_Mito)) +
   scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) + 
   geom_point(size = 1) + 
   theme_bw() + 
   theme(text = element_text(size = 20),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         axis.line = element_line(size = 1),
         axis.ticks = element_line(size = 1),
         legend.title = element_blank(),
         legend.key = element_blank())

ggplot(as.data.frame(C11_sce@reducedDims$PCA_coldata), 
       aes(x=PC1, y=PC2, color = C11_sce$log10_total_counts)) +
   scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) + 
   geom_point(size = 1) + 
   theme_bw() + 
   theme(text = element_text(size = 20),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         axis.line = element_line(size = 1),
         axis.ticks = element_line(size = 1),
         legend.title = element_blank(),
         legend.key = element_blank())


C11_sce$use <- (
   (C11_sce$pct_counts_Mito <= 10) &
      (C11_sce$log10_total_counts >= 3)
)
table(C11_sce$use)

ggplot(as.data.frame(C11_sce@reducedDims$PCA_coldata), 
       aes(x=PC1, y=PC2, color = as.factor(C11_sce$use))) +
   geom_point(size = 1) + 
   theme_bw() + 
   theme(text = element_text(size = 20),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         axis.line = element_line(size = 1),
         axis.ticks = element_line(size = 1),
         legend.title = element_blank(),
         legend.key = element_blank())

sce_C11_filtered = C11_sce[,C11_sce$use]
save(sce_C11_filtered, file = "C11_sce_filtered.RData")




# QC for C12
br.out <- barcodeRanks(counts(C12_sce))

plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")

set.seed(100)
e.out <- emptyDrops(counts(C12_sce))
e.out
is.cell <- e.out$FDR <= 0.05
sum(is.cell, na.rm=TRUE)

abline(h=br.out[br.out$rank == sum(is.cell, na.rm=TRUE)-1,]$total, col="red", lty=2)
abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen", "red"), 
       legend=c("knee", "inflection", "FDR_0.05"))

C12_sce <- C12_sce[,which(e.out$FDR <= 0.05)]
C12_sce@assays$data$logcounts <- log2(C12_sce@assays$data$counts + 1)

is.mito = rownames(C12_sce) %in% mtGenes[,1]
length(is.mito[is.mito == TRUE])
C12_sce = calculateQCMetrics(C12_sce, feature_controls=list(Mito=is.mito))


h <- hist(C12_sce$log10_total_counts, 
          breaks=100, col="grey80",
          xlab="Log-total UMI count")
cuts <- cut(h$breaks, c(0,3.0,Inf))
plot(h, col=c("red","white")[cuts])

h <- hist(C12_sce$pct_counts_Mito, 
          breaks=100, col="grey80",
          xlab="Proportion of reads in mitochondrial genes")
cuts <- cut(h$breaks, c(10,100,Inf))
plot(h, col=c("red","white")[cuts])


C12_sce <- runPCA(C12_sce, use_coldata=TRUE)
plotReducedDim(C12_sce, use_dimred = "PCA_coldata")

ggplot(as.data.frame(C12_sce@reducedDims$PCA_coldata), 
       aes(x=PC1, y=PC2, color = C12_sce$pct_counts_Mito)) +
   scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) + 
   geom_point(size = 1) + 
   theme_bw() + 
   theme(text = element_text(size = 20),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         axis.line = element_line(size = 1),
         axis.ticks = element_line(size = 1),
         legend.title = element_blank(),
         legend.key = element_blank())

ggplot(as.data.frame(C12_sce@reducedDims$PCA_coldata), 
       aes(x=PC1, y=PC2, color = C12_sce$log10_total_counts)) +
   scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) + 
   geom_point(size = 1) + 
   theme_bw() + 
   theme(text = element_text(size = 20),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         axis.line = element_line(size = 1),
         axis.ticks = element_line(size = 1),
         legend.title = element_blank(),
         legend.key = element_blank())


C12_sce$use <- (
   (C12_sce$pct_counts_Mito <= 10) &
      (C12_sce$log10_total_counts >= 3)
)
table(C12_sce$use)

ggplot(as.data.frame(C12_sce@reducedDims$PCA_coldata), 
       aes(x=PC1, y=PC2, color = as.factor(C12_sce$use))) +
   geom_point(size = 1) + 
   theme_bw() + 
   theme(text = element_text(size = 20),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         axis.line = element_line(size = 1),
         axis.ticks = element_line(size = 1),
         legend.title = element_blank(),
         legend.key = element_blank())

sce_C12_filtered = C12_sce[,C12_sce$use]
save(sce_C12_filtered, file = "C12_sce_filtered.RData")



var_list <- as.vector(ls())
sce_list <- grep("_filtered", var_list)

#assign gene names and cell barcodes into expression matrix
rownames(sce_C10_filtered@assays$data$counts) <- rownames(sce_C10_filtered)
colnames(sce_C10_filtered@assays$data$counts) <- paste0("PBS-", substr(sce_C10_filtered@colData$Barcode, 1, 16))
colnames(sce_C10_filtered) <- colnames(sce_C10_filtered@assays$data$counts)


rownames(sce_C11_filtered@assays$data$counts) <- rownames(sce_C11_filtered)
colnames(sce_C11_filtered@assays$data$counts) <- paste0("LPIMB19-", substr(sce_C11_filtered@colData$Barcode, 1, 16))
colnames(sce_C11_filtered) <- colnames(sce_C11_filtered@assays$data$counts)

rownames(sce_C12_filtered@assays$data$counts) <- rownames(sce_C12_filtered)
colnames(sce_C12_filtered@assays$data$counts) <- paste0("CPS-", substr(sce_C12_filtered@colData$Barcode, 1, 16))
colnames(sce_C12_filtered) <- colnames(sce_C12_filtered@assays$data$counts)


# Aggregate the matrices
for (i in sce_list){
  tmp <- get(var_list[i])
  mcols(tmp) <- mcols(C10_sce)
  assign(var_list[i], tmp)
  remove(tmp)
}


sce_filtered <- cbind(sce_C10_filtered, sce_C11_filtered)
sce_filtered <- cbind(sce_filtered, sce_C12_filtered)

sce_filtered <- sce_filtered[rowSums(as.matrix(counts(sce_filtered))) > 0,]
save(sce_filtered, file = "sce_filtered.RData")


# Normalization
clusters <- quickCluster(sce_filtered, method= "igraph")
sce_filtered <- computeSumFactors(sce_filtered, cluster=clusters)
sce_filtered <- normalize(sce_filtered)

sce_filtered_normalized <- sce_filtered
save(sce_filtered_normalized, file = "sce_filtered_normalized.RData")


#Determine highly variable genes
var.fit <- trendVar(sce_filtered, parametric=TRUE, use.spikes=F)
var.out <- decomposeVar(sce_filtered, var.fit)
hvg <- var.out[which(var.out$FDR <= 0.05),]
hvg <- hvg[order(hvg$bio, decreasing=TRUE),]

plot(var.out$mean, var.out$total, pch=16, cex=0.3, xlab="Mean log-expression",
     ylab="Variance of log-expression")
o <- order(var.out$mean)
lines(var.out$mean[o], var.out$tech[o], col="dodgerblue", lwd=2)
points(var.out$mean[var.out$FDR <= 0.05 & var.out$bio > 0.05],
       var.out$total[var.out$FDR <= 0.05 & var.out$bio > 0.05], pch=16, cex=0.3, col="red")


seurat <- as.Seurat(sce_filtered_normalized)
VariableFeatures(seurat) <- rownames(hvg)

#add condition on seurat
condition <- data.frame(rep(0,ncol(seurat)))
colnames(condition) <- "sample_condition"
condition[grep("PBS-", rownames(seurat@meta.data)),1] <- "PBS"
condition[grep("LPIMB19-", rownames(seurat@meta.data)),1] <- "LpIMB19"
condition[grep("CPS-", rownames(seurat@meta.data)),1] <- "RHP"
rownames(condition) <- colnames(seurat)

seurat@meta.data$sample_condition <- condition

#Dim reduction
seurat <- ScaleData(seurat, features = VariableFeatures(seurat))

seurat <- RunPCA(seurat, pcs.compute=50, weight.by.var = FALSE)
plot(seurat@reductions$pca@stdev)
PCA = 12

set.seed(19)
seurat <- FindNeighbors(seurat, reduction = "pca" ,dims = 1:PCA)
seurat <- FindClusters(seurat, dims = 1:PCA, save.SNN =T)
seurat <- RunUMAP(seurat, reduction = "pca", dims = 1:PCA, n.neighbors = 35L,
                        min.dist = 0.5, metric = "euclidean", seed.use = 42)
DimPlot(seurat, reduction = "umap", pt.size = 1, order = T)


#Remove non-cd8 t - expressing # Cd44+Il2rb+Klra1 / Pmel+Mlana
seurat <- seurat[, !(seurat@meta.data$seurat_clusters %in% c(5,10))]

seurat <- seurat[rowSums(seurat@assays$RNA@counts) != 0, ]
sce <- as.SingleCellExperiment(seurat)
