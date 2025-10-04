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


### True Counts
Truth <- read.csv("/home/isabel/Desktop/Thesis_publications_researchPapers/Thesis/Final_datasets_used/Codes_WD/Cell_Clustering/Synth_data/Truth.csv", row.names=NULL)
rownames(Truth) = Genes
colnames(Truth) = Cells
Truth = as.matrix(Truth)

#With scaling
heatmap(scale(Truth),col=gcol2,hclustfun=function(x)hclust(x,method="ward.D"))




## Noisy
Noisy_UMI_3000 <- read.csv("/home/isabel/Desktop/Thesis_publications_researchPapers/Thesis/Final_datasets_used/Codes_WD/Cell_Clustering/Synth_data/Noisy_UMI_3000.csv")
dim(Noisy_UMI_3000)
rownames(Noisy_UMI_3000) = Genes
colnames(Noisy_UMI_3000) = Cells
Noisy_UMI_3000 = as.matrix(Noisy_UMI_3000)
class(Noisy_UMI_3000)

#without scaling
heatmap(Noisy_UMI_3000,col=gcol2,hclustfun=function(x)hclust(x,method="ward.D"))
#With scaling
heatmap(scale(Noisy_UMI_3000),col=gcol2,hclustfun=function(x)hclust(x,method="ward.D"))



## SAVER
denoised_saver <- read.csv("~/Desktop/Thesis/Final_datasets_used/Codes_WD/Cell_Clustering/Synth_data/denoised_saver.csv")
View(denoised_saver)

rownames(denoised_saver)= Genes
colnames(denoised_saver) = Cells
class(denoised_saver)
denoised_saver=as.matrix(denoised_saver)

#With scaling
heatmap(scale(denoised_saver),col=gcol2,hclustfun=function(x)hclust(x,method="ward.D"))




## ENHANCE
denoised_enhance <- read.delim("/home/isabel/Desktop/Thesis_publications_researchPapers/Thesis/Final_datasets_used/Codes_WD/Cell_Clustering/Synth_data/denoised_enhance.tsv")
denoised_enhance = denoised_enhance[,-1]
rownames(denoised_enhance)=Genes
colnames(denoised_enhance)=Cells

denoised_enhance = as.matrix(denoised_enhance)

#With scaling
heatmap(scale(denoised_enhance),col=gcol2,hclustfun=function(x)hclust(x,method="ward.D"))




## MAGIC
denoised_magic <- read.csv("/home/isabel/Desktop/Thesis_publications_researchPapers/Thesis/Final_datasets_used/Codes_WD/Cell_Clustering/Synth_data/denoised_magic.csv")
View(denoised_magic)

denoised_magic = denoised_magic[,-1]
denoised_magic = t(denoised_magic)
denoised_magic = data.frame(denoised_magic)
rownames(denoised_magic)=Genes
colnames(denoised_magic)=Cells

dim(denoised_magic)
class(denoised_magic)
#denoised_magic = as.matrix(denoised_magic)

#With scaling
heatmap(scale(denoised_magic),col=gcol2,hclustfun=function(x)hclust(x,method="ward.D"))



## WAVDeSc
denoised_wavdesc <- read.csv("/home/isabel/Desktop/Thesis_publications_researchPapers/Thesis/Final_datasets_used/Codes_WD/Cell_Clustering/Synth_data/denoised_wavdesc_bior2.6D.csv", header=FALSE)
 View(denoised_wavdesc)


rownames(denoised_wavdesc)=Genes
colnames(denoised_wavdesc)=Cells

dim(denoised_wavdesc)
class(denoised_wavdesc)
denoised_wavdesc = as.matrix(denoised_wavdesc)

#With scaling
heatmap(scale(denoised_wavdesc),col=gcol2,hclustfun=function(x)hclust(x,method="ward.D"))

 ##### DCA
denoised_DCA <- read.delim("/home/isabel/Desktop/Thesis_publications_researchPapers/Thesis/Final_datasets_used/Codes_WD/Cell_Clustering/Synth_data/denoised_DCA.tsv", row.names=1)
View(denoised_DCA)

rownames(denoised_DCA)=Genes
colnames(denoised_DCA)=Cells

dim(denoised_DCA)
class(denoised_DCA)
denoised_DCA = as.matrix(denoised_DCA)

#With scaling
heatmap(scale(denoised_DCA),col=gcol2,hclustfun=function(x)hclust(x,method="ward.D"))






## WAVDeSc
Truth <- read.csv("/home/isabel/Desktop/Thesis_publications_researchPapers/Thesis/Final_datasets_used/Codes_WD/Cell_Clustering/Synth_data/Truth.csv")
   View(Truth)


rownames(Truth)=Genes
colnames(Truth)=Cells


Truth = as.matrix(Truth)

#With scaling
heatmap(scale(Truth),col=gcol2,hclustfun=function(x)hclust(x,method="ward.D"))



##### Putting all the heatmaps together


library(png)
library(jpeg)
library(grid)

Truth <- readPNG("/home/isabel/Desktop/Thesis_publications_researchPapers/Thesis/Final_datasets_used/Codes_WD/Cell_Clustering/cell_clustering_heatmap/Truth_heatmp.png")
Noisy <- readPNG("/home/isabel/Desktop/Thesis_publications_researchPapers/Thesis/Final_datasets_used/Codes_WD/Cell_Clustering/cell_clustering_heatmap/Noisy_heatmap.png")
Saver <- readPNG("/home/isabel/Desktop/Thesis_publications_researchPapers/Thesis/Final_datasets_used/Codes_WD/Cell_Clustering/cell_clustering_heatmap/SAVER_heatmp.png")
Magic <- readPNG("/home/isabel/Desktop/Thesis_publications_researchPapers/Thesis/Final_datasets_used/Codes_WD/Cell_Clustering/cell_clustering_heatmap/MAGIC_heatmp.png")
Enhance <- readPNG("/home/isabel/Desktop/Thesis_publications_researchPapers/Thesis/Final_datasets_used/Codes_WD/Cell_Clustering/cell_clustering_heatmap/ENHANCE_heatmp.png")
WAVDeSc <- readPNG("/home/isabel/Desktop/Thesis_publications_researchPapers/Thesis/Final_datasets_used/Codes_WD/Cell_Clustering/cell_clustering_heatmap/WAVDeSc.png")

# Display images
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 3)))  # 2 rows and 3 columns

# Plot the images
g1= rasterGrob(Truth, width = unit(1, "npc"), height = unit(1, "npc"))#,  vp = viewport(layout.pos.row = 1, layout.pos.col = 1)

g2=rasterGrob(Noisy, width = unit(1, "npc"), height = unit(1, "npc"))#,  vp = viewport(layout.pos.row = 1, layout.pos.col = 2)

g3 = rasterGrob(Saver, width = unit(1, "npc"), height = unit(1, "npc"))#,  vp = viewport(layout.pos.row = 1, layout.pos.col = 3)

g4 = rasterGrob(Magic, width = unit(1, "npc"), height = unit(1, "npc"))#, vp = viewport(layout.pos.row = 2, layout.pos.col = 1)

g5 = rasterGrob(Enhance, width = unit(1, "npc"), height = unit(1, "npc"))#, vp = viewport(layout.pos.row = 2, layout.pos.col = 2)

g6 = rasterGrob(WAVDeSc, width = unit(1, "npc"), height = unit(1, "npc"))#vp = viewport(layout.pos.row = 2, layout.pos.col = 3)

library(gridExtra)

grid.arrange(g1, g2, g3, g4, g5, g6, ncol = 3, 
             heights = unit(rep(3, 2), "cm"), widths = unit(rep(5, 3), "cm"))








library(cowplot)
library(magick)

# Arrange images in a grid with `cowplot`
plot_grid(
  ggdraw() + draw_image(Truth),
  ggdraw() + draw_image(Noisy),
  ggdraw() + draw_image(Saver),
  ggdraw() + draw_image(Enhance),
  ggdraw() + draw_image(Magic),
  ggdraw() + draw_image(WAVDeSc),
  ncol = 3
)











###### Real Data
reads <- read.delim("/home/isabel/Desktop/Thesis/Final_datasets_used/Real_data_nonUMI/tung/reads.csv")

original = counts_mat
geneNam =original[,1]

rownames(original) = geneNam
original =original[,-1]

original = as.matrix(original)
class(original)

#without scaling
heatmap(original,col=gcol2,hclustfun=function(x)hclust(x,method="ward.D"))
#With scaling
heatmap(scale(original),col=gcol2,hclustfun=function(x)hclust(x,method="ward.D"))




Real_saver <- read.csv("/home/isabel/Desktop/Thesis/Final_datasets_used/Real_data_nonUMI/Saver_nonUMI/denoised_Realsaver_nonUMI.csv")
rownames(Real_saver) = geneNam
Real_saver = as.matrix(Real_saver)

heatmap(scale(Real_saver),col=gcol2,hclustfun=function(x)hclust(x,method="ward.D"))


Real_saver <- read.csv("/home/isabel/Desktop/Thesis/Final_datasets_used/Real_data_nonUMI/Saver_nonUMI/denoised_Realsaver_nonUMI.csv")





