---
editor_options: 
  markdown: 
    wrap: 72
---

#calling require packages

library(SpatialExperiment)

library(STexampleData)

library(scater)

library(ggspavis)

#reading dataset

spe \<-
readRDS("C:/Users/sanaz/OneDrive/Documents/cosmx_lung_5_rep1_SPE.RDS")

#showing dimensions

dim(spe)

#View the names of the variables and the first few lines of data

head(colData(spe))

#Among the variables in this dataset, we select four variables, Width,
Height, Mean.PanCK, and Mean.DAPI, and in the next section, we calculate
the threshold value for them. With the following codes, the histogram of
all four variables can be drawn together:

par(mfrow=c(2,2))

hist(colData(spe)\$Width, breaks = 20)

hist(colData(spe)\$Height, breaks = 20)

hist(colData(spe)\$Mean.PanCK, breaks = 20)

hist(colData(spe)\$Mean.DAPI, breaks = 20)

![](Rplot04.png)

#Threshold value for the Width variable

#First, we redraw the histogram of this variable:

par(mfrow=c(1,1))

hist(colData(spe)\$Width, breaks = 20)

![](Rplot05.png)

#According to the horizontal axis, we consider the threshold value equal
to 200

#By choosing this value for the threshold, we draw the graph of the
Width variable in

#relation to the Area variable:

plotQC(spe, type = "scatter", metric_x = "Area", metric_y = "Width",
threshold_y = 200)![](Rplot.png)

#View the number of values that are larger and smaller than the
threshold value (200):

qc_Width = colData(spe)\$Width \> 200

table(qc_Width)

#Examining the pattern of discarded points:

colData(spe)\$qc_Width = qc_Width

plotQC(spe, type = "spots",in_tissue = NULL,discard = "qc_Width")

![](Rplot01.png)

#if consider the threshold 150:

qc_width_150 = colData(spe)\$Width \> 150

colData(spe)\$qc_width_150 = qc_width_150

plotQC(spe, type = "spots", in_tissue = NULL, discard = "qc_width_150")

![](Rplot02.png)

#Drawing layers:

plotSpots(spe, annotate = NULL,in_tissue = NULL, palette =
"libd_layer_colors")

![](Rplot03.png)

#Threshold value for Height variable

#First, we redraw the histogram of this variable:

hist(colData(spe)\$Height, breaks = 20)

![](Rplot06.png)

#According to the horizontal axis, we consider the threshold value equal
to 200. #By choosing this value for the threshold, #we draw the graph of
the Height variable in relation to the Area variable:

plotQC(spe, type = "scatter",

       metric_x = "Area", metric_y = "Height",

       threshold_y = 200)

![](Rplot07.png)

#View the number of values that are larger and smaller than the
threshold value (200):

qc_Height = colData(spe)\$Height \> 200

table(qc_Height)

#Examining the pattern of discarded points:

colData(spe)\$qc_Height = qc_Height

plotQC(spe, type = "spots",in_tissue = NULL,discard = "qc_Height")

![](Rplot08.png)

#If we mistakenly considered the threshold value to be 150, for example:

qc_Height_150 = colData(spe)\$Height \> 150

colData(spe)\$qc_Height_150 = qc_Height_150

plotQC(spe, type = "spots", in_tissue = NULL, discard = "qc_Height_150")

![](Rplot09.png)

#Drawing layers

plotSpots(spe, annotate = NULL,in_tissue = NULL, palette =
"libd_layer_colors")

![](Rplot10.png)

#Threshold value for Mean.PanCK variable

#First, we redraw the histogram of this variable:

hist(colData(spe)\$Mean.PanCK, breaks = 20)

![](Rplot11.png)

#According to the horizontal axis, we consider the threshold value equal
to 60000.

#By choosing this value for the threshold, we draw the graph of the

#Mean.PanCK variable against the Area variable:

plotQC(spe, type = "scatter",

       metric_x = "Area", metric_y = "Mean.PanCK",

       threshold_y = 60000)

![](Rplot12.png)

#View the number of values that are greater and less than the threshold
value (60000):

qc_Mean.PanCK = colData(spe)\$Mean.PanCK \> 60000

table(qc_Mean.PanCK)

#Examining the pattern of discarded points:

colData(spe)\$qc_Mean.PanCK = qc_Mean.PanCK

plotQC(spe, type = "spots",in_tissue = NULL,discard = "qc_Mean.PanCK")

![](Rplot13.png)

#If we consider the threshold value to be, for example, 50,000:

qc_Mean.PanCK_50000 = colData(spe)\$Mean.PanCK \> 50000

colData(spe)\$qc_Mean.PanCK_50000 = qc_Mean.PanCK_50000

plotQC(spe, type = "spots", in_tissue = NULL, discard =
"qc_Mean.PanCK_50000")

![](Rplot14.png)

#Drawing layers

plotSpots(spe, annotate = NULL,in_tissue = NULL, palette =
"libd_layer_colors")

![](Rplot15.png)

#Threshold value for Mean.DAPI variable

#First, we redraw the histogram of this variable:

hist(colData(spe)\$Mean.DAPI, breaks = 20)

![](Rplot16.png)

#According to the horizontal axis, we consider the threshold value equal
to 30000.

#By choosing this value for the threshold, #we draw the graph of the
Mean.DAPI variable against the Area variable:

plotQC(spe, type = "scatter",

       metric_x = "Area", metric_y = "Mean.DAPI",

       threshold_y = 30000)

![](Rplot17.png)

#View the number of values that are greater and less than the threshold
value (30000):

qc_Mean.DAPI = colData(spe)\$Mean.DAPI \> 30000

table(qc_Mean.DAPI)

#Examining the pattern of discarded points:

colData(spe)\$qc_Mean.DAPI = qc_Mean.DAPI

plotQC(spe, type = "spots",in_tissue = NULL,discard = "qc_Mean.DAPI")

![](Rplot18.png)

#If we consider the threshold value to be 20000 for example:

qc_Mean.DAPI_20000 = colData(spe)\$Mean.DAPI \> 20000

colData(spe)\$qc_Mean.DAPI_20000 = qc_Mean.DAPI_20000

plotQC(spe, type = "spots", in_tissue = NULL, discard =
"qc_Mean.DAPI_20000")

![](Rplot19.png)

#drawing layer

plotSpots(spe, annotate = NULL,in_tissue = NULL, palette =
"libd_layer_colors")

![](Rplot20.png)

#Removing points with poor quality

#Obtain the number of poor quality (outlier) points for each of the
variables discussed above:

apply(cbind(qc_Width, qc_Height, qc_Mean.PanCK, qc_Mean.DAPI), 2, sum)

#qc_Width qc_Height qc_Mean.PanCK qc_Mean.DAPI

\# 21 16 236 111

\# Combine poor quality points

discard = qc_Width \| qc_Height \| qc_Mean.PanCK \| qc_Mean.DAPI

table(discard)

#discard

\# FALSE TRUE

#99909 383

#We put the poor quality points on top of each other to draw their
graph:

colData(spe)\$discard = discard

#Draw a diagram of poor quality points:

plotQC(spe, type = "spots",in_tissue = NULL,discard = "discard")

![](Rplot21.png)

#Removing poor quality points and specifying new data dimensions after
removing these points:

spe = spe[, !colData(spe)\$discard]

dim(spe)

#980 99909

#houses with an area of zero and one

#Obtaining area abundances for some cells

tbl_cells_per_spot = table(colData(spe)\$Area)

tbl_cells_per_spot[1:13]

#4 20 24 28 32 36 40 44 48 52 56 60 64

\# 1 1 2 2 1 11 8 6 1 1 1 3 1

#Obtaining the area distribution for some cells

prop_cells_per_spot = round(tbl_cells_per_spot /
sum(tbl_cells_per_spot), 2)

prop_cells_per_spot[1:13]

#4 20 24 28 32 36 40 44 48 52 56 60 64

\# 0 0 0 0 0 0 0 0 0 0 0 0 0
