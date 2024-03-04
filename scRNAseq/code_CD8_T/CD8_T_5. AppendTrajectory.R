#add path
branch_probs <- read.csv("Palantir/palantir_res/palantir_branch_probs_til.csv", row.names = 1)
colnames(branch_probs) <- c("Exhausted", "Cytotoxic")
branch_probs <- branch_probs[colnames(seurat), ]
identical(rownames(branch_probs), colnames(seurat))


normWeights <- sweep(branch_probs, 1, FUN = "/",
                     STATS = apply(branch_probs, 1, sum))
head(normWeights)

set.seed(10)
wSamp <- apply(normWeights, 1, 
               function(prob){rmultinom(n = 1, prob = prob, size = 1)})
head(wSamp)

if(is.null(dim(wSamp))){
  wSamp <- matrix(wSamp, ncol = 1)
} else{
  wSamp <- t(wSamp)
}

colnames(wSamp) <- colnames(normWeights)
head(wSamp)

meta <- cbind(seurat@meta.data, wSamp)

seurat@meta.data <- meta

path_vec <- rep(NA, ncol(seurat))
names(path_vec) <- colnames(seurat)
path_vec[seurat$Cytotoxic == 1] <- "Cytotoxic"
path_vec[seurat$Exhausted == 1] <- "Exhausted"
seurat$path <- factor(path_vec, levels = c("Cytotoxic", "Exhausted"))
