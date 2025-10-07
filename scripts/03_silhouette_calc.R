# Calculate Silhouette Score
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
library(cluster)

silhouette_score <- function(Data, dataset_label){
  seurat_object = Seurat::CreateSeuratObject(counts = Data)
  seurat_object = NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
  all.genes <- rownames(seurat_object)
  seurat_object <- ScaleData(seurat_object, features = all.genes)
  seurat_object <- RunPCA(seurat_object, npcs = 50, features = all.genes, verbose = F)
  
  seurat_object <- FindNeighbors(seurat_object, dims = 1:20)
  seurat_object <- FindClusters(seurat_object, resolution = 0.5)
  predicted_labels <- Idents(seurat_object)
  seurat_object <- RunTSNE(seurat_object, dims = 1:10)
  p1 = DimPlot(seurat_object, reduction = "tsne", label = FALSE) 
  p1_customized <- p1 + ggtitle(dataset_label) + NoLegend() +
    theme(plot.title = element_text(color = "black", size = 20, hjust = 0.5, face = "bold"),
          plot.margin = unit(c(1,1,1,1), "lines"))
  
  # Calculate Silhouette Score
  silhouette_score <- silhouette(as.numeric(predicted_labels), dist(t(Data)))
  avg_silhouette_score <- mean(silhouette_score[, "sil_width"])
  
  return(list(plot = p1_customized, Silhouette_Score = avg_silhouette_score))
}


# NOISY
Noisy_UMI_3000 <- read.csv("Noisy_UMI_3000.csv")
results <- silhouette_score(Noisy_UMI_3000, "Noisy")
NSS =round(results$Silhouette_Score,3)

# MAGIC
denoised_magic <- read.csv("denoised_magic.csv", row.names = 1)
magic = t(denoised_magic)
resultsM<- silhouette_score(magic, "Magic")
MSS =round(resultsM$Silhouette_Score,3)


# SAVER
denoised_saver <- read_csv("denoised_saver.csv")
resultsS<- silhouette_score(denoised_saver, "SAVER")
SSS =round(resultsS$Silhouette_Score,3)


#ENHANCE
denoised_enhance <- read.table("denoised_enhance.tsv", row.names=1)
resultsE<- silhouette_score(denoised_enhance, "ENHANCE")
ESS =round(resultsE$Silhouette_Score,3)

#WAVDESC
denoised_wavdesc_bior2.6D <- read.csv("denoised_wavdesc_bior2.6D.csv", header=FALSE)
resultsW<- silhouette_score(denoised_wavdesc_bior2.6D, "WAVDESC")
WSS =round(resultsW$Silhouette_Score,3)

#DCA
denoised_DCA <- read.delim("denoised_DCA.tsv", row.names=1)
resultsD<- silhouette_score(denoised_DCA, "DCA")
DSS =round(resultsD$Silhouette_Score,3)

