#1	Quality Control
#Quality control (QC) procedures at the spot level aim to remove low-quality spots
#before further analysis. Low-quality spots can occur due to problems during library preparation
#or other experimental procedures. Examples include large proportions of dead cells due to cell
#damage during library preparation, and low mRNA capture efficiency due to inefficient reverse transcription or PCR amplification.
#These spots are usually removed prior to further analysis, since otherwise they tend to
#create problems during downstream analyses such as clustering. For example, problematic spots that are 
#not removed could show up as separate clusters, which may be misidentified as distinct cell types.
#Low-quality spots can be identified according to several characteristics, including:
#•	library size (i.e. total UMI counts per spot. It is denoted by “sum” in the data set.)
#•	number of expressed genes (i.e. number of genes with non-zero UMI counts per spot. It is denoted by “detected” in the data set.)
#Low library size or low number of expressed genes (or features) can indicate poor mRNA capture rates,
#e.g. due to cell damage and missing mRNAs, or low reaction efficiency.

#1-1	Load data
library(SpatialExperiment)
library(scater)
library(scran)
library(ggspavis)
spe <- readRDS("C:/Users/sanaz/OneDrive/Documents/cosmx_lung_5_rep1_SPE.RDS")
head(colData(spe))
##DataFrame with 6 rows and 19 columns

#1-2	Plot data
#As an initial check, plot the spatial coordinates (spots) 
#in x-y dimensions to check that the object has loaded correctly
#and that the orientation is as expected.
#Plot spatial coordinates (spots):
plotSpots(spe, in_tissue = NULL)


#1-3	Calculate QC metrics
#We calculate the QC metrics described above with a combination of 
#methods from the scater package. The QC metrics from scater can be
#calculated and added to the SpatialExperiment object as follows. 
#Calculate per-spot QC metrics (“sum” and “detected”) and store in colData:
spe = addPerCellQC(spe)
head(colData(spe))
##DataFrame with 6 rows and 22 columns


#Histogram histograms of QC metrics:
par(mfrow = c(1, 2))
hist(colData(spe)$sum, xlab = "sum", main = "UMIs per spot")
hist(colData(spe)$detected, xlab = "detected", main = "Genes per spot")

#1-4	Selecting thresholds
#The simplest option to apply the QC metrics is to select 
#thresholds for each metric, and remove any spots that do not meet the
#thresholds for one or more metrics. Exploratory visualizations can be 
#used to help select appropriate thresholds, which may differ depending on the dataset.
#Here, we use visualizations to select thresholds for several 
#QC metrics in our dataset: (i) library size, (ii) number of expressed genes (or features).

#1-4-1	Thresholds for library size (“sum”)
#Library size represents the total sum of UMI counts per spot. 
#This is included in the column labeled sum in the scater output.
#Histogram of library sizes:
par(mfrow=c(1,1))
hist(colData(spe)$sum, xlab = "sum",breaks=100, main = "UMIs per spot")
##The distribution is relatively smooth, but there is obvious issue such as a 
##spike at very low library sizes.
  
#Plot library size (“sum”) vs number of expressed genes (“detected”):
plotQC(spe, type = "scatter", 
       metric_x = "detected", metric_y = "sum", 
       threshold_y = 2)
##The horizontal line (argument threshold) shows our first
##guess at a possible filtering threshold for library size based on the histogram.
##The plot shows that setting a filtering threshold for library size
##(e.g. at the value shown) does not appear to select for any obvious
##biologically consistent group of spots.

#We set a relatively arbitrary threshold of 2 UMI counts per spot,
#and then check the number of spots below this threshold.
#Select QC threshold for library size:
qc_lib_size = colData(spe)$sum < 2
table(qc_lib_size)
##qc_lib_size
##FALSE   TRUE 
##100041    251 
colData(spe)$qc_lib_size = qc_lib_size
##Finally, we also check that the discarded
##spots do not have any obvious spatial pattern that correlates 
##with known biological features. Otherwise, removing these spots 
##could indicate that we have set the threshold too high, and are 
##removing biologically informative spots.

#Check spatial pattern of discarded spots:
plotQC(spe, type = "spots",in_tissue = NULL,discard = "qc_lib_size")
##discarded spots do not have any obvious spatial pattern


#1-4-2	Thresholds for Number of expressed genes (“detected”)
#The number of expressed genes (or features) refers to the number of
#genes with non-zero UMI counts per spot. This is stored in the column 
#detected in the scater output.
#We use a similar sequence of visualizations to choose a threshold for this QC metric.
#Histogram of numbers of expressed genes:
hist(colData(spe)$detected, xlab = "detected",breaks=100, main = "Genes per spot")
#Plot number of expressed genes (“detected”) vs library size (“sum”):
plotQC(spe, type = "scatter", 
       metric_x = "sum", metric_y = "detected", 
       threshold_y = 2)
##Based on the plots, we select a threshold of 2 for expressed genes per spot.

#Select QC threshold for number of expressed genes:
qc_detected = colData(spe)$detected < 2
table(qc_detected)
##qc_detected
##FALSE   TRUE 
##100033    259

#Check spatial pattern of discarded spots:
colData(spe)$qc_detected = qc_detected
plotQC(spe, type = "spots", in_tissue = NULL,discard = "qc_detected")

#1-4-3	Remove low-quality spots
#Now that we have calculated several QC metrics and
#selected thresholds for each one, we can combine the sets of 
#low-quality spots, and remove them from our object.
#We also check again that the combined set of discarded
#spots does not correspond to any obvious biologically relevant group of spots.
#Number of discarded spots for each metric:
apply(cbind(qc_lib_size, qc_detected), 2, sum)
##qc_lib_size qc_detected 
##251         259

#Combined set of discarded spots:
discard = qc_lib_size | qc_detected
table(discard)
##discard
##FALSE   TRUE 
##100033    259

#Store in object:
colData(spe)$discard = discard
#Check spatial pattern of combined set of discarded spots:
plotQC(spe, type = "spots",in_tissue = NULL,discard = "discard")

#Remove combined set of low-quality spots:
spe = spe[, !colData(spe)$discard]
dim(spe)
##[1]    980 100033


####2	Normalization
#2-1	Overview
#Here we apply normalization methods developed
#for scRNA-seq data, treating each spot as equivalent to one cell.
#2-2	Logcounts
#Calculate log-transformed normalized counts (abbreviated as “logcounts”) 
#using library size factors.
#We apply the methods implemented 
#in the scater and scran packages, which were originally developed for scRNA-seq
#data, making the assumption here that these methods can be applied to SRT data by 
#treating spots as equivalent to cells.
#We use the library size factors methodology since this is the simplest approach, 
#and can easily be applied to SRT data. Alternative approaches that are populare 
#for scRNA-seq data, including normalization by deconvolution, are more difficulty
#to justify in the context of spot-based SRT data since (i) 
#spots may contain multiple cells from more than one cell type, and (ii) datasets
#can contain multiple samples (e.g. multiple Visium slides, resulting
#in sample-specific clustering).
#Removing any spots with zero counts:
spe = spe[, colSums(counts(spe)) > 0]

#Calculate library size factors:
spe = computeLibraryFactors(spe)
summary(sizeFactors(spe))
hist(sizeFactors(spe), breaks = 100)
## Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
##0.006587 0.447949 0.820141 1.000000 1.350433 8.152006 

#Calculate logcounts and store in object:
spe = logNormCounts(spe)
#Check:
assayNames(spe)
##[1] "counts"    "logcounts"
dim(counts(spe))
##[1]    980 100033
dim(logcounts(spe))
##[1]    980 100033


#####3	Feature (gene) selection
#3-1	Overview
#Here we apply feature selection methods to
#identify highly variable genes (HVGs) or spatially variable genes (SVGs), 
#which can then be investigated individually or used as the input for further
#downstream analyses.

#3-2	Highly variable genes (HVGs)
#Scran methods are used to identify top highly variable genes (HVGs) 
#for defining major cell types. These methods, originally developed 
#for single-cell RNA sequencing data, do not consider spatial information.
#If the dataset mainly reflects spatial distributions of major cell types,
#HVGs may suffice, but if there are additional important spatial features,
#spatially variable genes (SVGs) may be more meaningful.

#Fit mean-variance relationship:
dec = modelGeneVar(spe)
fit = metadata(dec)

#visualize mean-variance relationship:
plot(fit$mean, fit$var, 
     xlab = "mean of log-expression", ylab = "variance of log-expression")
curve(fit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)

#Select top HVGs:
top_hvgs = getTopHVGs(dec, prop = 0.1)
length(top_hvgs)
##[1] 66

####4	Dimensionality reduction
#4-1	Overview
#In this section, we apply dimensionality reduction methods to visualize
#the data and to generate inputs for further downstream analyses.

#4-2	Principal component analysis (PCA)
#Principal component analysis (PCA) is applied to 
#top highly variable genes (HVGs) to reduce dataset dimensionality 
#and retain top 50 principal components for downstream analyses. 
#This reduces noise and improves computational efficiency. 
#The scater package provides a computationally efficient implementation.
#Compute PCA:
set.seed(123)
spe = runPCA(spe, subset_row = top_hvgs)
reducedDimNames(spe)
dim(reducedDim(spe, "PCA"))
##[1] "PCA"
##[1] 100033     50

#4-3	Uniform Manifold Approximation and Projection (UMAP)
#We also run UMAP on the set of top 50 PCs and retain the top 2 UMAP components, 
#which will be used for visualization purposes.
#Compute UMAP on top 50 PCs:
set.seed(123)
spe = runUMAP(spe, dimred = "PCA")
reducedDimNames(spe)
dim(reducedDim(spe, "UMAP"))
##[1] "PCA"  "UMAP"
##[1] 100033      2

#Update column names for easier plotting:
colnames(reducedDim(spe, "UMAP")) = paste0("UMAP", 1:2)



####4-4	Visualizations
#Create plots using ggspavis package's plotting functions,
#and in the next section on clustering, add cluster labels to reduced dimension plots.
#Plot top 2 PCA dimensions:
plotDimRed(spe, type = "PCA")

#Plot top 2 UMAP dimensions:
plotDimRed(spe, type = "UMAP")

##############5	Clustering
#5-1	Overview
#Clustering algorithms can be applied to ST data to 
#identify spatial domains, which are regions with consistent 
#gene expression profiles. These domains can be single cell types or
#a mixture of cell types. The optimal number of clusters depends on the
#biological context, and these domains can be further investigated using 
#differential expression testing to identify representative genes.

#5-2	Non-spatial clustering on HVGs
#The study uses standard clustering methods for single-cell RNA sequencing data, 
#focusing on molecular features like gene expression, to detect biologically 
#informative spatial distribution patterns of cell types in spatial data.

#Graph-based clustering:
set.seed(123)
k =10
g = buildSNNGraph(spe, k = k, use.dimred = "PCA")
g_walk = igraph::cluster_louvain(g)
clus = g_walk$membership
table(clus)
##clus
##1     2     3     4     5     6     7     8     9    10    11    12    13    14 
##12769  9132  4443  5652  6428  9850  2786 13151  7357  7609  2667  3230   421 14538 

#Store cluster labels in column 'label' in colData:
colLabels(spe) = factor(clus)

#Visualize tissue slide and PCA/UMAP clusters using ggspavis package functions.
#Results show clustering reproduces known biological structure, but not perfectly.
#Clusters are also not perfectly separated.
#Plot clusters in spatial x-y coordinates:
plotSpots(spe, annotate = colData(spe)$label, 
          palette = "libd_layer_colors",in_tissue = NULL)

#Plot clusters in PCA reduced dimensions:
plotDimRed(spe, type = "PCA", 
           annotate = colData(spe)$label, palette = "libd_layer_colors")

#Plot clusters in UMAP reduced dimensions:
plotDimRed(spe, type = "UMAP", 
           annotate = colData(spe)$label, palette = "libd_layer_colors")

#####6	Marker genes
#6-1	Overview
#In this section, we perform differential expression testing between
#clusters or spatial domains to identify representative marker genes for
#each cluster or spatial domain.
#6-2	Differential expression testing
#The findMarkers implementation in scran tests for 
#differential gene expression between clusters, using a binomial 
#test to select genes easier to interpret and validate experimentally.
#test for marker genes:
markers = findMarkers(spe, test = "binom", direction = "up")
#returns a list with one DataFrame per cluster:
markers
##List of length 14
##names(14): 1 2 3 4 5 6 7 8 9 10 11 12 13 14

#plot log-fold changes for one cluster over all other clusters (selecting cluster 1):
library(pheatmap)
interesting = markers[[1]]
best_set = interesting[interesting$Top <= 5, ]
logFCs = getMarkerEffects(best_set)
pheatmap(logFCs, breaks = seq(-5, 5, length.out = 101))

#plot log-transformed normalized expression of top genes for one cluster:
top_genes = head(rownames(interesting))
plotExpression(spe, x = "label", features = top_genes)






  