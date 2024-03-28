library(SingleCellExperiment)
library(DropletUtils)
library(scater)

sample_info_list <- list.files(dir)
plotdir <-dir
for( i in sample_info_list){
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
  sce
  assign(paste0(gsub("-","_",i),"_sce"),sce)
}
###Save Droplet QC SCE
for(i in sample_info_list){
  objectname = paste0(gsub("-","_",i),"_sce")
  sce = eval(parse(text=objectname))
  saveRDS(sce,paste0("data/",i,"sce_droplet100_QC.rds"))
}
dir.create("M:/dgcha_revision/plot/QC")
library(RColorBrewer)
source("D:/OneDrive - dgist.ac.kr/Function/QC_function.R")
####Set QC filter : log(counts)>2.5 mito <10
for(i in sample_info_list){
  objectname = paste0(gsub("-","_",i),"_sce")
  sce = eval(parse(text=objectname))
  sce1 <- qc_function(sce, objectname = objectname,
                      total_counts_cutoff = 3, mt_percent_cutoff = 10)
  sce.qc <- qc_function(sce, objectname = objectname,
                        total_counts_cutoff = 3, mt_percent_cutoff = 10, save = TRUE)
  assign(paste0(gsub("-","_",i),"_sce"), sce1)
  assign(paste0(gsub("-","_",i),"_sce", ".qc"), sce.qc)
}
colnames(E2_sce.qc) <- paste0(substring(colnames(E2_sce.qc),1,17),"PBS")
colnames(E4_sce.qc) <- paste0(substring(colnames(E4_sce.qc),1,17),"RHP")

saveRDS(E2_sce.qc,paste0("data/E2_umi3_mito10.rds"))
saveRDS(E4_sce.qc,paste0("data/E4_umi3_mito10.rds"))

