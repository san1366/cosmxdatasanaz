library(SpatialExperiment)
library(STexampleData)
library(scater)
library(scran)
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
spe <- computeLibraryFactors(spe)
summary(sizeFactors(spe))
hist(sizeFactors(spe), breaks = 20)
spe <- logNormCounts(spe)
assayNames(spe)
dim(counts(spe))
dim(logcounts(spe))
