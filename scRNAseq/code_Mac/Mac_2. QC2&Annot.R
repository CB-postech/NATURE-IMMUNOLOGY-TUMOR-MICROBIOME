.libPaths()
E2_sce.qc <- readRDS(paste0("data/E2_umi3_mito10.rds"))
E4_sce.qc <- readRDS(paste0("data/E4_umi3_mito10.rds"))
library(Seurat)
library(scater)
library(scran)
packageVersion("Seurat")
packageVersion("scater")
packageVersion("harmony")
packageVersion("harmony")
sce_merge <- cbind(E2_sce.qc,E4_sce.qc)
keep_feature <- rowSums(counts(sce_merge) != 0) > 0
sce_merge <- sce_merge[keep_feature, ]
### filtering genes
library(scran)
set.seed(123)
clusters <- quickCluster(sce_merge, method="igraph")
sce_merge <- computeSumFactors(sce_merge, cluster=clusters)
sce_merge <- logNormCounts(sce_merge)
dec <- modelGeneVar(sce_merge)
plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
curve(metadata(dec)$trend(x), col="blue", add=TRUE)
hvg <- getTopHVGs(dec,n=1000)
length(hvg)

library(Seurat)
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

Mac_seurat <- seurat_subset_recluster(seurat_bc, idents = c(4,9,13,14),fdr.threshold = 0.05,
                                      PCA = 20)
Mac_seurat_bc <-harmony_correction(Mac_seurat,batch = "sample",PCA = 20,assay = "originalexp")

qc_mac_seurat <- seurat_subset_recluster(Mac_seurat_bc, idents = c(7,9),invert = T,
                                         fdr.threshold = 0.05)
qc_mac_seurat_bc <-harmony_correction(qc_mac_seurat,batch = "sample",assay = "originalexp")

percent_celltype <- qc_mac_seurat_bc@meta.data %>% dplyr::group_by(sample,seurat_clusters) %>% dplyr::summarise(count=n()) %>% dplyr::mutate(perc=count/sum(count))
g1<- ggplot(percent_celltype, aes(x=sample, y=perc, fill=seurat_clusters)) +
  geom_bar(stat = "identity") + 
  # scale_fill_manual(values = brewer.pal(5,"Paired"))+
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


qc_mac_seurat_bc <- readRDS("paper_final/qc_mac_seurat_bc.rds")
DimPlot_eunseo(qc_mac_seurat_bc,label = T,group.by = "simple_celltype")
Lcn2_seurat <- seurat_subset_recluster(qc_mac_seurat_bc,cells = colnames(qc_mac_seurat_bc)[qc_mac_seurat_bc$simple_celltype=="Mac.Lcn2"],
                                       n=1000)
lcn2_cluster_marker <- FindAllMarkers(Lcn2_seurat,logfc.threshold = 0.4)
Lcn2_seurat_bc <- harmony_correction(Lcn2_seurat,batch = "sample")
lcn2_cluster_marker_bc <- FindAllMarkers(Lcn2_seurat_bc,logfc.threshold = 0.4)

Lcn2_seurat_bc$celltype <- "Inflammatory macrophage"
Lcn2_seurat_bc$celltype[Lcn2_seurat_bc$seurat_clusters %in% c(4,5,0,3,6)] <- "Regulatory TAM"
Lcn2_seurat_bc <- ScaleData(Lcn2_seurat_bc,features = rownames(Lcn2_seurat_bc))

qc_mac_seurat_bc$celltype <- "Monocytic macrophage"
qc_mac_seurat_bc$celltype[colnames(Lcn2_seurat_bc)] <- as.character(Lcn2_seurat_bc$celltype)

library(RColorBrewer)
percent_celltype <- qc_mac_seurat_bc@meta.data %>% dplyr::group_by(sample,celltype) %>% dplyr::summarise(count=n()) %>% dplyr::mutate(perc=count/sum(count))
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
g1