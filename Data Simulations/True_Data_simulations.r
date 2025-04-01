## Load packages
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


## Default simulation parameters from paper

## 1. Simulating Truecount (500 by 300)
phyla1 <- Phyla5()-5 cell populations
ngenes <- 500
true_counts_res <- SimulateTrueCounts(ncells_total=300, min_popsize=50, i_minpop=2, ngenes=ngenes, nevf=10,
                                      evf_type="discrete",n_de_evf=9, vary="s", Sigma=0.4, phyla=Phyla5(), randseed=0)
Truth1 = true_counts_res[["counts"]]
View(Truth1)



## 2. Simulating Truecount (6000by3000)- 3 cell populations
phyla2 <- read.tree(system.file("extdata", "Newick_AB.txt", package = "SymSim"))
plot(phyla2)

#phyla2 <- read.tree(system.file("extdata", "Newick_ABCDE.txt", package = "SymSim"))
ngenes <- 6000
true_counts_res <- SimulateTrueCounts(ncells_total=3000, min_popsize=800, i_minpop=3, ngenes=ngenes, nevf=20,
                                      evf_type="discrete",n_de_evf=9, vary="s", gene_effect_prob = 0.1, gene_effects_sd = 0.5, Sigma=0.5, phyla=phyla2, randseed=10)
Truth2 = true_counts_res[["counts"]]
View(Truth2)




## 3. Simulating Truecount (10000by3000)- 5 cell populations
phyla1 <- Phyla5()
phyla2 <- read.tree(system.file("extdata", "Newick_ABCDE.txt", package = "SymSim"))
ngenes <- 10000
true_counts_res <- SimulateTrueCounts(ncells_total=3000, min_popsize=600, ngenes=ngenes, nevf=60,
                                      evf_type="discrete",n_de_evf=9, vary="s", gene_effect_prob = 0.1, gene_effects_sd = 0.5, Sigma=0.5, phyla=Phyla5(), randseed=10)
Truth1 = true_counts_res[["counts"]]





## plotting clusters for True data
true_counts_res_dis <- true_counts_res
tsne_true_counts <- PlotTsne(meta=true_counts_res[[3]], data=log2(true_counts_res[[1]]+1), evf_type="discrete", 
                             n_pc=20,label='pop', saving = F, plotname="discrete populations (true counts)")
tsne_true_counts[[2]]

## plotting clusters for True data using seurate
tmp_so <- summarizedToSeurat(true_counts_res)
clustering_result = calculate_clustering(tmp_so)


