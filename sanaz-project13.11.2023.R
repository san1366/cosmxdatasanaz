#load packages
library(SpatialExperiment)
library(STexampleData)
library(scater)
library(ggspavis)
#load data
spe <- readRDS("C:/Users/sanaz/OneDrive/Documents/cosmx_lung_5_rep1_SPE.RDS")
dim(spe)
head(colData(spe))
#Hisogram
par(mfrow=c(2,2))
hist(colData(spe)$Width, breaks = 20)
hist(colData(spe)$Height, breaks = 20)
hist(colData(spe)$Mean.PanCK, breaks = 20)
hist(colData(spe)$Mean.DAPI, breaks = 20)
#	Selecting thresholds:
#	Thresholds for Width
#Histogram for width
par(mfrow=c(1,1))
hist(colData(spe)$Width, breaks = 20)
plotQC(spe, type = "scatter", 
       metric_x = "Area", metric_y = "Width", 
       threshold_y = 200)
qc_Width = colData(spe)$Width > 200
table(qc_Width)
colData(spe)$qc_Width = qc_Width
plotQC(spe, type = "spots",in_tissue = NULL,discard = "qc_Width")
qc_width_150 = colData(spe)$Width > 150
colData(spe)$qc_width_150 = qc_width_150
plotQC(spe, type = "spots", in_tissue = NULL, discard = "qc_width_150")
#drawing layers
plotSpots(spe, annotate = NULL,in_tissue = NULL, palette = "libd_layer_colors")

#	Thresholds for Height
#Histogram
hist(colData(spe)$Height, breaks = 20)

plotQC(spe, type = "scatter", 
       metric_x = "Area", metric_y = "Height", 
       threshold_y = 200)

qc_Height = colData(spe)$Height > 200
table(qc_Height)

colData(spe)$qc_Height = qc_Height
plotQC(spe, type = "spots",in_tissue = NULL,discard = "qc_Height")


qc_Height_150 = colData(spe)$Height > 150
colData(spe)$qc_Height_150 = qc_Height_150
plotQC(spe, type = "spots", in_tissue = NULL, discard = "qc_Height_150")

plotSpots(spe, annotate = NULL,in_tissue = NULL, palette = "libd_layer_colors")

#	Thresholds for Mean.PanCK
#histogram
hist(colData(spe)$Mean.PanCK, breaks = 20)

plotQC(spe, type = "scatter", 
       metric_x = "Area", metric_y = "Mean.PanCK", 
       threshold_y = 60000)

qc_Mean.PanCK = colData(spe)$Mean.PanCK > 60000
table(qc_Mean.PanCK)

colData(spe)$qc_Mean.PanCK = qc_Mean.PanCK
plotQC(spe, type = "spots",in_tissue = NULL,discard = "qc_Mean.PanCK")

qc_Mean.PanCK_50000 = colData(spe)$Mean.PanCK > 50000
colData(spe)$qc_Mean.PanCK_50000 = qc_Mean.PanCK_50000
plotQC(spe, type = "spots", in_tissue = NULL, discard = "qc_Mean.PanCK_50000")

plotSpots(spe, annotate = NULL,in_tissue = NULL, palette = "libd_layer_colors")

#	Thresholds for Mean.DAPI
#histogram
hist(colData(spe)$Mean.DAPI, breaks = 20)


plotQC(spe, type = "scatter", 
       metric_x = "Area", metric_y = "Mean.DAPI", 
       threshold_y = 30000)

qc_Mean.DAPI = colData(spe)$Mean.DAPI > 30000
table(qc_Mean.DAPI)


colData(spe)$qc_Mean.DAPI = qc_Mean.DAPI
plotQC(spe, type = "spots",in_tissue = NULL,discard = "qc_Mean.DAPI")


qc_Mean.DAPI_20000 = colData(spe)$Mean.DAPI > 20000
colData(spe)$qc_Mean.DAPI_20000 = qc_Mean.DAPI_20000
plotQC(spe, type = "spots", in_tissue = NULL, discard = "qc_Mean.DAPI_20000")

plotSpots(spe, annotate = NULL,in_tissue = NULL, palette = "libd_layer_colors")


#	Remove low-quality spots
apply(cbind(qc_Width, qc_Height, qc_Mean.PanCK, qc_Mean.DAPI), 2, sum)

discard = qc_Width | qc_Height | qc_Mean.PanCK | qc_Mean.DAPI
table(discard)


colData(spe)$discard = discard
plotQC(spe, type = "spots",in_tissue = NULL,discard = "discard")

spe = spe[, !colData(spe)$discard]
dim(spe)

#	Zero-cell and single-cell spots
tbl_cells_per_spot = table(colData(spe)$Area)
tbl_cells_per_spot[1:13]

prop_cells_per_spot = round(tbl_cells_per_spot / sum(tbl_cells_per_spot), 2)
prop_cells_per_spot[1:13]


#	Normalization
#Removing any spots with non-positive counts:
spe = spe[, colSums(counts(spe)) > 0]
#Calculate library size factors:
spe <- computeLibraryFactors(spe)
summary(sizeFactors(spe))
hist(sizeFactors(spe), breaks = 20)

#Calculate logcounts and store in object:
spe <- logNormCounts(spe)
#Check:
assayNames(spe)
dim(counts(spe))
dim(logcounts(spe))

#Feature selection
#Fit mean-variance relationship:
library(scran)
dec = modelGeneVar(spe)
fit = metadata(dec)
#visualize mean-variance relationship:
plot(fit$mean, fit$var, 
     xlab = "mean of log-expression", ylab = "variance of log-expression")
curve(fit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)
#Select top HVGs:
top_hvgs = getTopHVGs(dec, prop = 0.1)
length(top_hvgs)
  


#	Dimensionality reduction
#Compute PCA:
set.seed(123)
spe = runPCA(spe, subset_row = top_hvgs)
reducedDimNames(spe)
dim(reducedDim(spe, "PCA"))

#Compute UMAP on top 50 PCs:
set.seed(123)
spe = runUMAP(spe, dimred = "PCA")
reducedDimNames(spe)
dim(reducedDim(spe, "UMAP"))
#Update column names for easier plotting:
colnames(reducedDim(spe, "UMAP")) = paste0("UMAP", 1:2)
#Plot top 2 PCA dimensions:
plotDimRed(spe, type = "PCA")
#Plot top 2 UMAP dimensions:
plotDimRed(spe, type = "UMAP")

#	Clustering
#spe1 is a Subset of data for clustering:
spe1=spe[1:200,1:10000]

#Graph-based clustering:
set.seed(123)
k =10
g = buildSNNGraph(spe1, k = k, use.dimred = "PCA")
g_walk = igraph::cluster_walktrap(g)
clus = g_walk$membership
table(clus)
#Store cluster labels in column 'label' in colData:
colLabels(spe1) = factor(clus)
#Plot clusters in spatial x-y coordinates:
plotSpots(spe1, annotate = colData(spe1)$label, 
          palette = "libd_layer_colors",in_tissue = NULL)
#Plot clusters in PCA reduced dimensions:
plotDimRed(spe1, type = "PCA", 
           annotate = NULL, palette = "libd_layer_colors")
#Plot clusters in UMAP reduced dimensions:
spe1=spe[1:200,1:10000]
colLabels(spe1) = factor(clus)
plotDimRed(spe1, type = "UMAP", 
           annotate = colData(spe1)$label, palette = "libd_layer_colors")



