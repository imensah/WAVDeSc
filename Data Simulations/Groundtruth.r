### We simulated scRNA-seq data with different dimensions of rows and columns
#using the symsim software. Firstly we simulated the groundtruth to represent
#the truth data and then visualised the results using clusters.
#We then used then added on noise to the simulated groundtruth to produce the
#noisy data.

## Load all packages needed for the analysis
library(SymSim)
library(stats)
library(coop)
library(ggplot2)
library(cowplot)
library(fdm2id)
library(Metrics)
library(ggpubr)
library(Seurat)
library(SeuratObject)
library(SummarizedExperiment)


## NB: Here simulations are done with some of the the default simulation
#parameters used in the symsim paper.

## 1. Simulating Truecount with dimensions 500 genes by 300 cells with
# five(5) cell populations
phyla1 <- Phyla5()
ngenes <- 500
true_counts_res <- SimulateTrueCounts(
  ncells_total = 300,
  min_popsize = 50,
  i_minpop = 2,
  ngenes = ngenes,
  nevf = 10,
  evf_type = "discrete",
  n_de_evf = 9,
  vary = "s",
  Sigma = 0.4,
  phyla = Phyla5(),
  randseed = 0
)

#Extraction and checking the true count matrix
groundtruth_1 <- true_counts_res[["counts"]]
View(groundtruth_1)



## 2. Simulating Truecount with dimensions 6000 genes by 3000 cells
# with three(3) cell populations
phyla2 <- read.tree(system.file("extdata", "Newick_AB.txt", package = "SymSim"))
plot(phyla2)

ngenes <- 6000
true_counts_res <- SimulateTrueCounts(
  ncells_total = 3000,
  min_popsize = 800,
  i_minpop = 3,
  ngenes = ngenes,
  nevf = 20,
  evf_type = "discrete",
  n_de_evf = 9,
  vary = "s",
  gene_effect_prob = 0.1,
  gene_effects_sd = 0.5,
  Sigma = 0.5,
  phyla = phyla2,
  randseed = 10
)

#Extraction and checking the true count matrix for ground truth two
groundruth_2 <- true_counts_res[["counts"]]
View(groundruth_2)




## 3. Simulating Truecount with dimensions 10000 genes by 3000 cells
# with five cell populations
phyla1 <- Phyla5()

phyla2 <- read.tree(
  system.file("extdata", "Newick_ABCDE.txt", package = "SymSim")
)

ngenes <- 10000

true_counts_res <- SimulateTrueCounts(
  ncells_total = 3000,
  min_popsize = 600,
  ngenes = ngenes,
  nevf = 60,
  evf_type = "discrete",
  n_de_evf = 9,
  vary = "s",
  gene_effect_prob = 0.1,
  gene_effects_sd = 0.5,
  Sigma = 0.5,
  phyla = Phyla5(),
  randseed = 10
)


groundtruth_3 <- true_counts_res[["counts"]]
View(groundruth_3)





## plotting clusters for True data
true_counts_res_dis <- true_counts_res

tsne_true_counts <- PlotTsne(
  meta = true_counts_res[[3]],
  data = log2(true_counts_res[[1]] + 1),
  evf_type = "discrete",
  n_pc = 20,
  label = 'pop',
  saving = FALSE,
  plotname = "discrete populations (true counts)"
)

tsne_true_counts[[2]]

## plotting clusters for True data using seurate
tmp_so <- summarizedToSeurat(true_counts_res)
clustering_result <- calculate_clustering(tmp_so)