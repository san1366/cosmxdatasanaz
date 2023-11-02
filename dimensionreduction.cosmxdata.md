\#**Code to run steps from the previous chapters to generate the SpatialExperiment object required for this chapter:**

library(SpatialExperiment)

library(STexampleData)

library(scater)

library(scran)

library(ggspavis)

spe \<- readRDS("C:/Users/sanaz/OneDrive/Documents/cosmx_lung_5_rep1_SPE.RDS")

qc_Width = colData(spe)\$Width \> 200

colData(spe)\$qc_Width = qc_Width

qc_Height = colData(spe)\$Height \> 200

colData(spe)\$qc_Height = qc_Height

qc_Mean.PanCK = colData(spe)\$Mean.PanCK \> 60000

colData(spe)\$qc_Mean.PanCK = qc_Mean.PanCK

qc_Mean.DAPI = colData(spe)\$Mean.DAPI \> 30000

colData(spe)\$qc_Mean.DAPI = qc_Mean.DAPI

apply(cbind(qc_Width, qc_Height, qc_Mean.PanCK, qc_Mean.DAPI), 2, sum)

discard = qc_Width \| qc_Height \| qc_Mean.PanCK \| qc_Mean.DAPI

table(discard)

colData(spe)\$discard = discard

spe = spe[, !colData(spe)\$discard]

spe = spe[, colSums(counts(spe)) \> 0]

spe = computeLibraryFactors(spe)

spe = logNormCounts(spe)

dec = modelGeneVar(spe)

top_hvgs = getTopHVGs(dec, prop = 0.1)

\#**Compute PCA:**

set.seed(123)

spe = runPCA(spe, subset_row = top_hvgs)

reducedDimNames(spe)

dim(reducedDim(spe, "PCA"))

\#''PCA''

#99772 50

#Compute UMAP on top 50 PCs:

set.seed(123)

spe = runUMAP(spe, dimred = "PCA")

reducedDimNames(spe)

dim(reducedDim(spe, "UMAP"))

\#''PCA'' ''UMAP"

#99772 2

#Update column names for easier plotting:

colnames(reducedDim(spe, "UMAP")) = paste0("UMAP", 1:2)

\#**Plot top 2 PCA dimensions:**

plotDimRed(spe, type = "PCA")

![]()

![](dimension red1.png)

![]()

**#Plot top 2 UMAP dimensions:**

plotDimRed(spe, type = "UMAP")

![](dimension red2.png)
