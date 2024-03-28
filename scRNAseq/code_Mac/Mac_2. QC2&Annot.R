library(Seurat)
library(scater)
library(scran)
library(dplyr)
library(RColorBrewer)

harmony_correction <- function(seurat_obj, batch=NULL, PCA=15, assay="RNA"){
  library(Seurat)
  library(harmony)
  library(dplyr)

  PCA <- PCA
  set.seed(123456)
  seurat_H <- seurat_obj %>% 
    RunHarmony(batch, plot_convergence = F, max.iter.harmony = 100,assay.use =assay )
  seurat_H <- seurat_H %>% 
    RunUMAP(reduction = "harmony", dims = 1:PCA, seed.use = 42) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:PCA) %>% 
    FindClusters(resolution = 0.8) %>% 
    identity()
  return(seurat_H)
}

seurat_subset_recluster <- function(seurat_obj, cells=NULL, idents=NULL, invert=FALSE,
                                    PCA=15,n = NULL, var.threshold = 0, fdr.threshold = NULL){
  library(Seurat)
  library(scran)
  library(scater)

  seurat_subset <- subset(seurat_obj, cells=cells, idents=idents, invert=invert)
  sce <- as.SingleCellExperiment(seurat_subset)

  set.seed(123)
  clusters <- quickCluster(sce, method="igraph")
  table(clusters)
  sce <- computeSumFactors(sce, cluster=clusters)
  sce <- logNormCounts(sce)

  dec <- modelGeneVar(sce)
  hvg <- getTopHVGs(dec, n=n,
                    var.threshold = var.threshold, fdr.field = "FDR",
                    fdr.threshold = fdr.threshold)
  
  print(paste0("The number of hvg is ",length(hvg)))
  
  seurat <- as.Seurat(sce)
  VariableFeatures(seurat) <- hvg
  seurat <- ScaleData(seurat)
  
  seurat@reductions$PCA_coldata <- NULL
  seurat@reductions$pca <- NULL
  seurat@reductions$umap <- NULL
  seurat@reductions$PCA <- NULL
  seurat@reductions$UMAP <- NULL
  
  PCA = PCA
  
  seurat<- RunPCA(seurat, npcs = PCA)
  seurat<- FindNeighbors(seurat, reduction="pca", dims= 1:PCA)
  seurat<- FindClusters(seurat)
  seurat<- RunUMAP(seurat, dims = 1:PCA, seed.use = 42)
  
  return(seurat)
}


E2_sce.qc <- readRDS(paste0("data/E2_umi3_mito10.rds"))
E4_sce.qc <- readRDS(paste0("data/E4_umi3_mito10.rds"))

sce_merge <- cbind(E2_sce.qc,E4_sce.qc)
keep_feature <- rowSums(counts(sce_merge) != 0) > 0
sce_merge <- sce_merge[keep_feature, ]


### filtering genes
set.seed(123)
clusters <- quickCluster(sce_merge, method="igraph")
sce_merge <- computeSumFactors(sce_merge, cluster=clusters)
sce_merge <- logNormCounts(sce_merge)
dec <- modelGeneVar(sce_merge)
plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
curve(metadata(dec)$trend(x), col="blue", add=TRUE)
hvg <- getTopHVGs(dec,n=1000)


seurat <- as.Seurat(sce_merge)
VariableFeatures(seurat) <- hvg
seurat <- ScaleData(seurat)
seurat@reductions$PCA_coldata <- NULL

PCA = 20
seurat<- RunPCA(seurat, npcs = PCA)
seurat<- FindNeighbors(seurat, reduction="pca", dims= 1:PCA)
seurat<- FindClusters(seurat)
seurat<- RunUMAP(seurat, dims = 1:PCA, seed.use = 42)


seurat$sample <- substring(colnames(seurat),18)
seurat_bc <-harmony_correction(seurat,batch = "sample",PCA = 20,assay = "originalexp")


Mac_seurat <- seurat_subset_recluster(seurat_bc, idents = c(4,9,13,14), fdr.threshold = 0.05, PCA = 20)
Mac_seurat_bc <-harmony_correction(Mac_seurat, batch = "sample", PCA = 20, assay = "originalexp")

qc_mac_seurat <- seurat_subset_recluster(Mac_seurat_bc, idents = c(7,9), invert = T, fdr.threshold = 0.05)
qc_mac_seurat_bc <-harmony_correction(qc_mac_seurat, batch = "sample", assay = "originalexp")

Lcn2_seurat <- seurat_subset_recluster(qc_mac_seurat_bc, cells = colnames(qc_mac_seurat_bc)[qc_mac_seurat_bc$simple_celltype == "Mac.Lcn2"], n=1000)
Lcn2_seurat_bc <- harmony_correction(Lcn2_seurat,batch = "sample")

Lcn2_seurat_bc$celltype <- "Inflammatory macrophage"
Lcn2_seurat_bc$celltype[Lcn2_seurat_bc$seurat_clusters %in% c(4,5,0,3,6)] <- "Regulatory TAM"
Lcn2_seurat_bc <- ScaleData(Lcn2_seurat_bc, features = rownames(Lcn2_seurat_bc))

qc_mac_seurat_bc$celltype <- "Monocytic macrophage"
qc_mac_seurat_bc$celltype[colnames(Lcn2_seurat_bc)] <- as.character(Lcn2_seurat_bc$celltype)


percent_celltype <- qc_mac_seurat_bc@meta.data %>%
                        group_by(sample,celltype) %>%
                        summarise(count=n()) %>%
                        mutate(perc=count/sum(count))


g1<- ggplot(percent_celltype, aes(x=sample, y=perc*100, fill=celltype)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = brewer.pal(4,"Paired"))+
  theme(  # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    text = element_text(size=10),
    axis.text.x = element_text(size = 10,colour = "black",angle = 90))