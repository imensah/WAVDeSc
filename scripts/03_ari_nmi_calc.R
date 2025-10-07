# --FOR SIMULATED DATA
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


############ ARI and NMI calculation
ari_nmi <- function(Data, dataset_label, true_labels) {
  seurat_object <- Seurat::CreateSeuratObject(counts = Data)
  seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
  all.genes <- rownames(seurat_object)
  seurat_object <- ScaleData(seurat_object, features = all.genes)
  seurat_object <- RunPCA(seurat_object, npcs = 50, features = all.genes, verbose = FALSE)
  
  seurat_object <- FindNeighbors(seurat_object, dims = 1:20)
  seurat_object <- FindClusters(seurat_object, resolution = 0.5)
  predicted_labels <- as.numeric(Idents(seurat_object))  # Ensure this is a numeric vector
  
  seurat_object <- RunTSNE(seurat_object, dims = 1:10)
  p1 <- DimPlot(seurat_object, reduction = "tsne", label = FALSE)
  p1_customized <- p1 + ggtitle(dataset_label) + NoLegend() +
    theme(plot.title = element_text(color = "black", size = 20, hjust = 0.5, face = "bold"),
          plot.margin = unit(c(1,1,1,1), "lines"))
  
  # Calculate ARI and NMI
  if (length(true_labels) == length(predicted_labels)) {
    ari_score <- mclust::adjustedRandIndex(true_labels, predicted_labels)
    nmi_score <- aricode::NMI(true_labels, predicted_labels)
  } else {
    ari_score <- NA  # Set to NA if lengths differ
    nmi_score <- NA  # Set to NA if lengths differ
  }
  
  # Return plot and metrics
  return(list(plot = p1_customized, predicted_labels, ARI = ari_score, NMI = nmi_score))
}


# Assuming true_labels is loaded from a CSV and needs adjustment
true_labels <- read.csv("true_labels.csv", stringsAsFactors = FALSE)
true_labels_vector <- as.numeric(true_labels$pop)  

# NOISY
Noisy <- read.csv("Noisy_UMI_3000.csv")
NOISY <- ari_nmi(Noisy, "Noisy", true_labels_vector)
AR1 = NOISY$ARI
round(AR1, 3)
NMI =NOISY$NMI
round(NMI,3)


## MAGIC
denoised_magic <- read.csv("denoised_magic.csv", row.names=1)
denoised_magic = t(denoised_magic)
dim(denoised_magic)
MAGIC <- ari_nmi(denoised_magic, "MAGIC", true_labels_vector)
AR1 = MAGIC$ARI
round(AR1, 3)
NMI =MAGIC$NMI
round(NMI,3)

### SAVER
denoised_saver <- read.csv("denoised_saver.csv")
SAVER <- ari_nmi(denoised_saver, "SAVER", true_labels_vector)
AR1 = SAVER$ARI
round(AR1, 3)
NMI =SAVER$NMI
round(NMI,3)


### ENHANCE
denoised_enhance <- read.table("denoised_enhance.tsv", row.names=1)
Enhance <- ari_nmi(denoised_enhance, "ENHANCE", true_labels_vector)
AR1 = Enhance$ARI
round(AR1, 3)
NMI =Enhance$NMI
round(NMI,3)


## WAVDESC
denoised_wavdesc_bior2.6D <- read.csv("denoised_wavdesc_bior2.6D.csv", header=FALSE)
WAVDESC <- ari_nmi(denoised_wavdesc_bior2.6D, "WAVDESC", true_labels_vector)
AR1 = WAVDESC$ARI
round(AR1, 1)
NMI =WAVDESC$NMI
round(NMI,1)

## DCA
denoised_DCA <- read.table("denoised_DCA.tsv", row.names=1)
DCA <- ari_nmi(denoised_DCA, "DCA", true_labels_vector)
AR1 = DCA$ARI
round(AR1, 1)
NMI =DCA$NMI
round(NMI,1)


# --FOR PBMC DATA (REAL DATA)
Cell_clustering <- function(Data, dataset_label, true_labels) {
  seurat_object <- Seurat::CreateSeuratObject(counts = Data)
  seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
  all.genes <- rownames(seurat_object)
  seurat_object <- ScaleData(seurat_object, features = all.genes)
  seurat_object <- RunPCA(seurat_object, npcs = 50, features = all.genes, verbose = FALSE)
  
  seurat_object <- FindNeighbors(seurat_object, dims = 1:20)
  seurat_object <- FindClusters(seurat_object, resolution = 0.5)
  predicted_labels <- as.numeric(Idents(seurat_object))  # Ensure this is a numeric vector
  
  #seurat_object <- RunTSNE(seurat_object, dims = 1:10)
  seurat_object <- RunTSNE(seurat_object, dims = 1:10,perplexity = 30, check_duplicates = FALSE)
  p1 <- DimPlot(seurat_object, reduction = "tsne", label = FALSE)
  p1_customized <- p1 + ggtitle(dataset_label) + NoLegend() +
    theme(plot.title = element_text(color = "black", size = 20, hjust = 0.5, face = "bold"),
          plot.margin = unit(c(1,1,1,1), "lines"))
  
  # Calculate ARI and NMI
  if (length(true_labels) == length(predicted_labels)) {
    ari_score <- mclust::adjustedRandIndex(true_labels, predicted_labels)
    nmi_score <- aricode::NMI(true_labels, predicted_labels)
  } else {
    ari_score <- NA  # Set to NA if lengths differ
    nmi_score <- NA  # Set to NA if lengths differ
  }
  
  # Calculate Silhouette Score
  # silhouette_score <- silhouette(as.numeric(predicted_labels), dist(t(Data)))
  # avg_silhouette_score <- mean(silhouette_score[, "sil_width"])
  
  
  # Return plot and metrics
  return(list(plot = p1_customized, predicted_labels, ARI = ari_score, NMI = nmi_score))#,Silhouette_Score = avg_silhouette_score))
}

pbmc.4k_expression <- read.delim("pbmc-4k_expression.tsv", row.names=1)

# Assuming true_labels is loaded from a CSV and needs adjustment
True_labels_pbmc <- read.csv("True_labels_pbmc.csv", row.names=1)
true_labels_vector <- as.numeric(True_labels_pbmc$Cluster)  

# NOISY
results <- Cell_clustering(pbmc.4k_expression, "Original", true_labels_vector)
AR1 = results$ARI
round(AR1, 3)
NMI =results$NMI
round(NMI,3)
#SS =results$Silhouette_Score
#round(SS,3)


## MAGIC
MAGIC_denoised_pbmc <- read.csv("MAGIC_denoised_pbmc.csv", row.names=1)
denoised_magic = t(MAGIC_denoised_pbmc)
dim(denoised_magic)

MAGIC <- Cell_clustering(denoised_magic, "MAGIC", true_labels_vector)
AR1 = MAGIC$ARI
round(AR1, 3)
NMI =MAGIC$NMI
round(NMI,3)

### SAVER
SAVER_denoised_pbmc <- read.csv("SAVER_denoised_pbmc.csv", row.names=1)
SAVER <- Cell_clustering(SAVER_denoised_pbmc, "SAVER", true_labels_vector)
AR1 = SAVER$ARI
round(AR1, 3)
NMI =SAVER$NMI
round(NMI,3)


### ENHANCE
ENHANCE_denoised_pbmc <- read.delim("ENHANCE_denoised_pbmc.tsv", row.names=1)
Enhance <- Cell_clustering(ENHANCE_denoised_pbmc, "ENHANCE", true_labels_vector)
AR1 = Enhance$ARI
round(AR1, 3)
NMI =Enhance$NMI
round(NMI,3)



## DCA
DCA_denoised_pbmc <- read.delim("DCA_denoised_pbmc.tsv", row.names=1)
DCA <- Cell_clustering(DCA_denoised_pbmc, "DCA", true_labels_vector)
AR1 = DCA$ARI
round(AR1,3)
NMI =DCA$NMI
round(NMI,3)



## WAVDESC
wavdesc_denoised_pbmc <- read.csv("wavdesc_denoised_pbmc.csv", header=FALSE)
WAVDESC <- Cell_clustering(wavdesc_denoised_pbmc, "WAVDESC", true_labels_vector)
AR1 = WAVDESC$ARI
round(AR1, 3)
NMI =WAVDESC$NMI
round(NMI,3)

