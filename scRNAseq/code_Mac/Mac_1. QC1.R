library(SingleCellExperiment)
library(DropletUtils)
library(scater)
library(RColorBrewer)


qc_function <- function(sce, total_counts_cutoff=0, total_features_cutoff=0,
                        mt_percent_cutoff=0, save = FALSE, objectname=NULL){
  if(sum(is.null(objectname)) == 1){
    objectname = deparse(substitute(sce))
    objectname = gsub("raw", "", objectname)
  }
  
  path = paste0(plotdir, "/QC/", strsplit(objectname, "_")[[1]][1], "_qc/")
  if(sum(dir.exists(path)) == 0){
    dir.create(path)
  }
  path = paste0(path, objectname, "/")
  if(sum(dir.exists(path)) == 0){
    dir.create(path)
  }

  load("ensemblGenes2022-02-21.RData")
  mtGenes = ensemblGenes[ensemblGenes[,3] == "MT",]
  is.mito = rownames(sce) %in% mtGenes[,2]
  print(paste0("There is ", length(is.mito[is.mito == TRUE]), " mitochondrial genes"))
  
  per.cell <- perCellQCMetrics(
    sce,
    subsets = list(MT=is.mito),
    percent_top = c(50, 100, 200, 500)
  )
  
  colData(sce) <- cbind(colData(sce), per.cell)
  
  sce$log10_sum = log10(sce$sum + 1)
  sce$log10_detected = log10(sce$detected + 1)
  
  ### plot function
  coldata = c("log10_sum", "log10_detected", "subsets_MT_percent")
  columnName = c("log10 total umi counts", "log10 total feature counts", "Mitochondrial percentage (%)")
  cutoffs = c(total_counts_cutoff, total_features_cutoff, mt_percent_cutoff)
  
  # run pca
  vector <-  c(unique(colnames(sce@colData)))[-c(1,2)]
  sce <- runColDataPCA(sce, ncomponents=5, variables = vector)
  
  # hist
  for(i in 1:length(coldata)){
    qc_hist_function(sce, v=cutoffs[i], coldata[i], columnName[i], path=path)
  }
  
  df <- as.data.frame(sce@int_colData$reducedDims$PCA_coldata)
  df$log10_sum <- sce$log10_sum
  df$log10_detected <- sce$log10_detected
  df$subsets_MT_percent <- sce$subsets_MT_percent
  
  # pca plot
  for(i in 1:length(coldata)){
    qc_dimplot_function(df, coldata[i], columnName[i], path=path)
  }
  
  filter_by_total_counts = sce$log10_sum > total_counts_cutoff
  filter_by_feature_counts = sce$log10_detected > total_features_cutoff
  filter_by_mt_percent = sce$subsets_MT_percent < mt_percent_cutoff
  
  sce$use <- (
    filter_by_feature_counts &
      filter_by_mt_percent &
      filter_by_total_counts
  )
  print(table(sce$use))
  
  ### save function
  if(sum(save) == 1){
    
    sce.qc <- sce[, sce$use]
    
    g <-plotReducedDim(sce, "PCA_coldata",
                   colour_by = "use",
                   size_by = "detected") +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.line=element_blank(),
            axis.ticks=element_blank(),
            axis.title = element_blank(),
            panel.border = element_blank(),
            legend.text=element_text(size=13),
            legend.key=element_blank(),
            axis.text = element_blank()) 
    ggsave(paste0(path, "dimplot_use.png"), width=6, height=5,g)

    print(paste0("sce.qc is generated. There are ", ncol(sce.qc), " cells finally."))

    return(sce.qc)
  } else{
    return(sce)
  }
}


sample_info_list <- list.files(dir)
plotdir <- dir
for(i in sample_info_list){
  dir <- paste0(dir,i,"/raw_feature_bc_matrix")
  sce <- read10xCounts(dir)
  rownames(sce) = uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)
  colnames(sce) = sce$Barcode
  
  my.counts=counts(sce)
  br.out <- barcodeRanks(my.counts)
  png(paste0(i,"_100.png"))
  plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
  o <- order(br.out$rank)
  lines(br.out$rank[o], br.out$fitted[o], col="red")
  abline(h=100, col="red", lty=2)
  abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
  abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
  legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
         legend=c("knee", "inflection"))
  dev.off()

  set.seed(100)
  e.out <- emptyDrops(my.counts,lower = 100)
  is.cell <- e.out$FDR <= 0.05
  print(paste0("The number of QC positive cells of ",i," is ",sum(is.cell, na.rm=TRUE)))
  is.cell[is.na(is.cell)] <- FALSE
  
  names(is.cell) <- colnames(sce)
  sce$cells_kept <- is.cell
  
  sce <- sce[, sce$cells_kept == T]

  assign(paste0(gsub("-","_",i),"_sce"), sce)
}


###Save Droplet QC SCE
for(i in sample_info_list){
  objectname = paste0(gsub("-","_",i),"_sce")
  sce = eval(parse(text=objectname))
  saveRDS(sce,paste0("data/",i,"sce_droplet100_QC.rds"))
}

####Set QC filter : log(counts)>2.5 mito <10
for(i in sample_info_list){
  objectname = paste0(gsub("-","_",i),"_sce")
  sce = eval(parse(text=objectname))

  sce1 <- qc_function(sce, objectname = objectname, total_counts_cutoff = 3, mt_percent_cutoff = 10)
  sce.qc <- qc_function(sce, objectname = objectname, total_counts_cutoff = 3, mt_percent_cutoff = 10, save = TRUE)

  assign(paste0(gsub("-","_",i),"_sce"), sce1)
  assign(paste0(gsub("-","_",i),"_sce", ".qc"), sce.qc)
}

colnames(E2_sce.qc) <- paste0(substring(colnames(E2_sce.qc),1,17),"PBS")
colnames(E4_sce.qc) <- paste0(substring(colnames(E4_sce.qc),1,17),"RHP")

saveRDS(E2_sce.qc, "data/E2_umi3_mito10.rds")
saveRDS(E4_sce.qc, "data/E4_umi3_mito10.rds")