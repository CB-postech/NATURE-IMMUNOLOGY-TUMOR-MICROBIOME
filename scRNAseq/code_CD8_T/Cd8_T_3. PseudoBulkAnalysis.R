library(DESeq2)

pseudobulk_seurat <- SplitObject(seurat, split.by = "CellType")
for(i in names(pseudobulk_seurat)){
  pseudobulk_seurat[[i]] <- SplitObject(pseudobulk_seurat[[i]], split.by = "sample_condition")
}

#Raw counts
for(i in names(pseudobulk_seurat)){
  for(j in names(pseudobulk_seurat[[i]])){
    pseudobulk_seurat[[i]][[j]] <- pseudobulk_seurat[[i]][[j]]@assays$RNA@counts
    tmp_gene_list <- rownames(pseudobulk_seurat[[i]][[j]])
    pseudobulk_seurat[[i]][[j]] <- rowSums(pseudobulk_seurat[[i]][[j]])
    names(pseudobulk_seurat[[i]][[j]]) <- tmp_gene_list
  }
}

#Check gene order
for(i in names(pseudobulk_seurat)){
  for(j in names(pseudobulk_seurat[[i]])){
    print(identical(names(pseudobulk_seurat[[i]][[j]]), tmp_gene_list))
  }
}

df_bulk = data.frame(row.names = tmp_gene_list)

#Add sample identities
for(i in c(1:length(names(pseudobulk_seurat)))){
  for(j in c(1:length(names(pseudobulk_seurat[[i]])))){
    df_bulk[[paste0(names(pseudobulk_seurat[[i]])[j], "_", names(pseudobulk_seurat)[i])]] <- pseudobulk_seurat[[i]][[j]]
  }
}

mtx_bulk <- as.matrix(df_bulk)
celltype_vec <- factor(rep(levels(seurat$CellType), 3), levels = levels(seurat$CellType))
mtx_bulk <- mtx_bulk[,c(7, 8, 9, 10, 11, 12, 4, 5, 6, 13, 14, 15, 1, 2, 3)]
annot <- data.frame(sample_condition = factor(rep(c("PBS", "LpIMB19", "RHP"), 5), levels = c("PBS", "LpIMB19", "RHP")),
                    CellType = sort(celltype_vec)
)
dds <- DESeqDataSetFromMatrix(mtx_bulk, colData = annot, design = ~ sample_condition)

#Normalization
dds <- normTransform(dds)
head(dds@assays@data@listData[[1]])

pcaData <- plotPCA(dds, intgroup = c("sample_condition", "CellType"), returnData = T)
percentVar <- round(100*attr(pcaData, "percentVar"))
p <- ggplot(pcaData, aes(PC1, PC2, color = sample_condition, shape = CellType)) +
  geom_point(data = pcaData, aes(PC1, PC2, shape = CellType, col = sample_condition, fill = sample_condition), size = 3) +
  scale_color_manual(values=c("#FF0000", "#000000", "#8C8C8C"), labels = c("PBS", "LpIMB19", "RHP")) +
  scale_fill_manual(values=c("#FF0000", "#000000", "#8C8C8C"), labels = c("PBS", "LpIMB19", "RHP")) +
  scale_shape_manual(labels = c("Naive", "Central Memory", "Effector", "Effector Memory", "Exhausted"),
                     values=c(21, 22, 23, 24, 25)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.text = element_text(colour = "#000000"),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())
