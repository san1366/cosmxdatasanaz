\#**Normalization**

\#**Code to run steps from the previous chapters to generate the SpatialExperiment object required for this chapter:**

library(SpatialExperiment)

library(STexampleData)

library(scater)

library(scran)

spe = readRDS("cosmx_lung_5_rep1_SPE.RDS")

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

\#**Removing any spots with non-positive counts:**

spe = spe[, colSums(counts(spe)) \> 0]

\#**Calculate library size factors:**

spe \<- computeLibraryFactors(spe)

summary(sizeFactors(spe))

hist(sizeFactors(spe), breaks = 20)

#Min. 1st Qu. Median Mean 3rd Qu. Max.

\# 0.003301 0.445653 0.821983 1.000000 1.350164 8.170311

![](nor,1.png)

#Calculate logcounts and store in object:

spe \<- logNormCounts(spe)

\#**Check:**

assayNames(spe)

dim(counts(spe))

dim(logcounts(spe))

\#"count" "logcounts"

#980 99772

#980 99772
