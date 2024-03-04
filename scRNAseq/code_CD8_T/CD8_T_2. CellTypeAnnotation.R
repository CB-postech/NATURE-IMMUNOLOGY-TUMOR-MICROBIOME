Cluster2CellType <- function(SeuratObj, ClusterSlot = "seurat_clusters", # name of slot as input
                             CellTypeSlot = "CellType", # name of cell type slot to be returned
                             QueryList = list(), # list of (name of cell type) to (name of clusters)
                             ReturnPlot = T){
  
  celltype_vec <- c(rep(NA, ncol(SeuratObj)))
  names(celltype_vec) <- colnames(SeuratObj)
  
  for(i in names(QueryList)){
    for(j in QueryList[[i]]){
      celltype_vec[SeuratObj[[ClusterSlot]] == j] <- i
    }
  }
  
  SeuratObj[[CellTypeSlot]] <- factor(celltype_vec, levels = names(QueryList))
  
  if(ncol(SeuratObj) == sum(table(celltype_vec))){
    print(table(celltype_vec))
    cat(paste0("Total #Cells: ", sum(table(celltype_vec)), "\n"))
  }else{
    stop(paste0("Total #Cells Differ: ", ncol(SeuratObj), ", ", sum(table(celltype_vec)), "\n"))
  }

  remove(celltype_vec)
  
  if(ReturnPlot == T){
    p <- DimPlot(SeuratObj, group.by = CellTypeSlot, label = T, label.size = 7, pt.size = 1) +
          theme(axis.line=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                panel.background=element_blank(),
                panel.border=element_blank(),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                plot.background=element_blank()
      )
    print(p)
  }
  return(SeuratObj)
}

#Normalization
clusters <- quickCluster(sce, method = "igraph")
sce <- computeSumFactors(sce, clusters = clusters)
sce <- logNormCounts(sce)

#Determine HVGs
dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, fdr.threshold = 0.05)

#Dim reduction
seurat <- as.Seurat(sce)
seurat@reductions$PCA <- NULL
seurat@reductions$UMAP <- NULL
VariableFeatures(seurat) <- hvg
seurat <- ScaleData(seurat, features = hvg)
seurat <- RunPCA(seurat, weight.by.var = F, npcs = 50)
plot(seurat@reductions$pca@stdev)
PCA = 10

set.seed(10)
seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:PCA)
seurat <- FindClusters(seurat, resolution = 0.8)
seurat <- RunUMAP(seurat, reduction = "pca", dims = 1:PCA, n.neighbors = 20)


query_list <- list("Naive" = c(1, 2, 6), "Central Memory" = c(3, 4, 5), "Effector" = c(0), "Effector Memory" = c(9), "Exhausted" = c(7, 8))
seurat <- Cluster2CellType(seurat, QueryList = query_list, ReturnPlot = T)


#Add Module Score for Cytotoxicity geneset
geneset_act <- c("Prf1", "Ifng", "Gnly", "Nkg7", "Gzmb", "Gzma", "Gzmh", "Klrk1", "Klrb1", "Klrd1", "Ctsw", "Cst7")
geneset_act <- geneset_act[geneset_act %in% rownames(seurat)]
geneset_act <- list(geneset_act)

seurat <- AddModuleScore(seurat, features = geneset_act, name = "geneset_act")
