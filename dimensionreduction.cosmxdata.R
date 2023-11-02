library(SpatialExperiment)
library(STexampleData)
library(scater)
library(scran)
library(ggspavis)
spe <- readRDS("C:/Users/sanaz/OneDrive/Documents/cosmx_lung_5_rep1_SPE.RDS")
qc_Width = colData(spe)$Width > 200
colData(spe)$qc_Width = qc_Width
qc_Height = colData(spe)$Height > 200
colData(spe)$qc_Height = qc_Height
qc_Mean.PanCK = colData(spe)$Mean.PanCK > 60000
colData(spe)$qc_Mean.PanCK = qc_Mean.PanCK
qc_Mean.DAPI = colData(spe)$Mean.DAPI > 30000
colData(spe)$qc_Mean.DAPI = qc_Mean.DAPI
apply(cbind(qc_Width, qc_Height, qc_Mean.PanCK, qc_Mean.DAPI), 2, sum)
discard = qc_Width | qc_Height | qc_Mean.PanCK | qc_Mean.DAPI
table(discard)
colData(spe)$discard = discard
spe = spe[, !colData(spe)$discard]
spe = spe[, colSums(counts(spe)) > 0]
spe = computeLibraryFactors(spe)
spe = logNormCounts(spe)
dec = modelGeneVar(spe)
top_hvgs = getTopHVGs(dec, prop = 0.1)
set.seed(123)
spe = runPCA(spe, subset_row = top_hvgs)
reducedDimNames(spe)
dim(reducedDim(spe, "PCA"))
set.seed(123)
spe = runUMAP(spe, dimred = "PCA")
reducedDimNames(spe)
dim(reducedDim(spe, "UMAP"))
colnames(reducedDim(spe, "UMAP")) = paste0("UMAP", 1:2)
plotDimRed(spe, type = "PCA")

plotDimRed(spe, type = "UMAP")




