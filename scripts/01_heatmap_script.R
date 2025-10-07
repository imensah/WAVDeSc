library(ggplot2)
library(kknn)
library(GGally)
library(umap)
library(Rtsne)
library(igraph)
library(huge)

Genes = paste("Genes",1:2000, sep = " ")
Cells = paste("Cells", 1:3000, sep = " ")

#cluster heatmap - biclustering
aa = grep("grey",colors())
bb = grep("green",colors())
cc = grep("red",colors())
gcol2 = colors()[c(aa[1:2],bb[1:25],cc[1:50])]
#gcol2 = colors()[c(aa,bb,cc)]


## plotting heatmaps ##

## -- Noisy
Noisy_UMI_3000 <- read.csv("/home/isabel/Desktop/Dataset4Heatmap/Noisy_UMI_3000.csv")
dim(Noisy_UMI_3000)
rownames(Noisy_UMI_3000) = Genes
colnames(Noisy_UMI_3000) = Cells
Noisy_UMI_3000 = as.matrix(Noisy_UMI_3000)
class(Noisy_UMI_3000)

# Heatmap With scaling
heatmap(scale(Noisy_UMI_3000),col=gcol2,hclustfun=function(x)hclust(x,method="ward.D"))



## -- SAVER
denoised_saver <- read.csv("/home/isabel/Desktop/Dataset4Heatmap/denoised_saver.csv")
View(denoised_saver)

rownames(denoised_saver)= Genes
colnames(denoised_saver) = Cells
class(denoised_saver)
denoised_saver=as.matrix(denoised_saver)

# Heatmap With scaling
heatmap(scale(denoised_saver),col=gcol2,hclustfun=function(x)hclust(x,method="ward.D"))



## -- ENHANCE
denoised_enhance <- read.delim("/home/isabel/Desktop/Dataset4Heatmap/denoised_enhance.tsv")
denoised_enhance = denoised_enhance[,-1]
rownames(denoised_enhance)=Genes
colnames(denoised_enhance)=Cells
denoised_enhance = as.matrix(denoised_enhance)

#With scaling
heatmap(scale(denoised_enhance),col=gcol2,hclustfun=function(x)hclust(x,method="ward.D"))



## MAGIC
denoised_magic <- read.csv("/home/isabel/Desktop/Dataset4Heatmap/denoised_magic.csv")
View(denoised_magic)
denoised_magic = denoised_magic[,-1]
denoised_magic = t(denoised_magic)
denoised_magic = data.frame(denoised_magic)
rownames(denoised_magic)=Genes
colnames(denoised_magic)=Cells

dim(denoised_magic)
class(denoised_magic)
denoised_magic = as.matrix(denoised_magic)

#With scaling
heatmap(scale(denoised_magic),col=gcol2,hclustfun=function(x)hclust(x,method="ward.D"))



## --WAVDeSc
denoised_wavdesc <- read.csv("/home/isabel/Desktop/Dataset4Heatmap/denoised_wavdesc_bior2.6D.csv", header=FALSE)
View(denoised_wavdesc)
rownames(denoised_wavdesc)=Genes
colnames(denoised_wavdesc)=Cells
dim(denoised_wavdesc)
class(denoised_wavdesc)
denoised_wavdesc = as.matrix(denoised_wavdesc)

#With scaling
heatmap(scale(denoised_wavdesc),col=gcol2,hclustfun=function(x)hclust(x,method="ward.D"))



## --DCA
denoised_DCA <- read.delim("/home/isabel/Desktop/Dataset4Heatmap/denoised_DCA.tsv", row.names=1)
View(denoised_DCA)
rownames(denoised_DCA)=Genes
colnames(denoised_DCA)=Cells
dim(denoised_DCA)
class(denoised_DCA)
denoised_DCA = as.matrix(denoised_DCA)

# With scaling
heatmap(scale(denoised_DCA),col=gcol2,hclustfun=function(x)hclust(x,method="ward.D"))


