# -- FOR SIMULATED DATA
library(pcaMethods)
library(SC3)
library(scater)
library(SingleCellExperiment)
library(pheatmap)
library(mclust)
library(igraph)
library(scran)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(readr)
library(aricode)


setwd("../clustering")
getwd()
list.files()

################ clustering for 2000genes by 3000cells #######################

#Write a function for clustering
Cell_clustering = function(Data, dataset_label){
  seurat_object = Seurat::CreateSeuratObject(counts = Data)
  seurat_object = NormalizeData(seurat_object,  normalization.method = "LogNormalize", 
                                scale.factor = 10000)
  all.genes <- rownames(seurat_object)
  seurat_object <- ScaleData(seurat_object, features = all.genes)
  seurat_object <- RunPCA(seurat_object,npcs = 50, features = all.genes, verbose = F)
  
  # Clustering
  seurat_object <- FindNeighbors(seurat_object, dims = 1:20)
  seurat_object <- FindClusters(seurat_object, resolution = 0.5)
  seurat_object <- RunTSNE(seurat_object, dims = 1:10)
  
  p1 = DimPlot(seurat_object, reduction = "tsne", label = FALSE)
  
  p1_customized <- p1 +  
    ggtitle(dataset_label) +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      axis.title = element_text(size = 14)
    )
  
  return(p1_customized) 
}

# -- noisy
Noisy = read.csv(file= "Noisy_UMI_3000.csv")
P1 = Cell_clustering(Noisy,"NOISY (3000 cells, UMI)")

Cell_clustering = function(Data, dataset_label){
  seurat_object = Seurat::CreateSeuratObject(counts = Data)
  seurat_object = NormalizeData(seurat_object,  normalization.method = "LogNormalize", 
                                scale.factor = 10000)
  all.genes <- rownames(seurat_object)
  seurat_object <- ScaleData(seurat_object, features = all.genes)
  seurat_object <- RunPCA(seurat_object,npcs = 50, features = all.genes, verbose = F)
  
  # Clustering
  seurat_object <- FindNeighbors(seurat_object, dims = 1:20)
  seurat_object <- FindClusters(seurat_object, resolution = 0.5)
  seurat_object <- RunTSNE(seurat_object, dims = 1:10)
  
  p1 = DimPlot(seurat_object, reduction = "tsne", label = FALSE) + NoLegend()
  
  p1_customized <- p1 +  
    ggtitle(dataset_label) +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      axis.title = element_text(size = 14)
    )
  
  return(p1_customized) 
}

# -- magic
Magic = read.csv(file= "denoised_magic.csv",row.names = 1)
Magic = t(Magic)
dim(Magic)
P2 = Cell_clustering(Magic,"MAGIC")

# --saver
Saver = read.csv(file= "denoised_saver.csv")
P3 = Cell_clustering(Saver,"SAVER")

# --enhance
Enhance = read.table(file= "denoised_enhance.tsv",row.names = 1)
P4 = Cell_clustering(Enhance,"ENHANCE")

# --dca
denoised_DCA <- read.table("denoised_DCA.tsv", row.names=1)
P5= Cell_clustering(denoised_DCA,"DCA")

# --wavdesc
Wavdesc = read.csv(file= "denoised_wavdesc_bior2.6D.csv", header = F)
P6 = Cell_clustering(Wavdesc,"WAVDESC")


combined_plot <- (P1 | P2 | P3) /
                 (P4 | P5 | P6)




# -- FOR PBMC DATA (REAL DATA)
################ clustering for PBMC #######################

setwd("../clustering/")
getwd()
list.files()

# Write a function for clustering
Cell_clustering = function(Data,dataset_label){
  seurat_object = Seurat::CreateSeuratObject(counts = Data)
  seurat_object = NormalizeData(seurat_object,  normalization.method = "LogNormalize", 
                                scale.factor = 10000)
  all.genes <- rownames(seurat_object)
  seurat_object <- ScaleData(seurat_object, features = all.genes)
  seurat_object <- RunPCA(seurat_object,npcs = 50, features = all.genes, verbose = F)
  #print(seurat_object[["pca"]], dims = 1:5, nfeatures = 5)
  
  
  #Clustering
  
  seurat_object <- FindNeighbors(seurat_object, dims = 1:20)
  seurat_object <- FindClusters(seurat_object, resolution = 0.5)
  #head(Idents(seurat_object), 5)
  # seurat_object <- RunTSNE(seurat_object, dims = 1:10)
  #In instances of duplication use
  seurat_object <- RunTSNE(seurat_object, dims = 1:10,perplexity = 30, check_duplicates = FALSE)
  p1=DimPlot(seurat_object, reduction = "tsne",label = FALSE) + NoLegend()
  p1_customized <- p1 +  
    ggtitle(dataset_label) +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      axis.title = element_text(size = 14)
    )
  
  return(p1_customized) 
}

# --Noisy
Noisy <- read.delim(file = "pbmc-4k_expression.tsv", header = TRUE, sep = "\t")
pbmcCol =  Noisy[,1]
rownames(Noisy) =pbmcCol
Noisy =Noisy[,-1]
P1 = Cell_clustering(Noisy,"PBMC 4K (4334 cells, 10x)")


# --Magic
Magic = read.csv(file= "MAGIC_denoised_pbmc.csv",row.names = 1)
Magic = t(Magic)
P2 = Cell_clustering(Magic,"MAGIC")


# --Saver
Saver = read.csv(file= "SAVER_denoised_pbmc.csv")
pbmcSav =  Saver[,1]
rownames(Saver) =pbmcSav
Saver =Saver[,-1]
Saver = as.matrix(Saver)
P3 = Cell_clustering(Saver,"SAVER")

# --Enhance
ENHANCE_denoised_pbmc <- read.csv("ENHANCE_denoised_pbmc.tsv", row.names = 1)
P4 <- Cell_clustering(ENHANCE_denoised_pbmc, "ENHANCE")

# --DCA
DCA_denoised_pbmc <- read.delim("DCA_denoised_pbmc.tsv", row.names=1)
P5 = Cell_clustering(DCA_denoised_pbmc,"DCA")

# --WAVDESC
wavdesc_denoised_pbmc <- read.csv("wavdesc_denoised_pbmc.csv", header=FALSE)
P6 = Cell_clustering(wavdesc_denoised_pbmc,"WAVDESC")


combined_plot <- (P1 | P2 | P3) /
  (P4 | P5 | P6)

