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


##1a Simulating Noisycounts using default parameters
set.seed(111)
data(gene_len_pool)
gene_len <- sample(gene_len_pool, ngenes, replace = FALSE)

observed_counts <- True2ObservedCounts(true_counts=true_counts_res[[1]], meta_cell=true_counts_res[[3]], 
                                       protocol="nonUMI", alpha_mean=0.1, alpha_sd=0.05, gene_len=gene_len,
                                       depth_mean=1e5, depth_sd=3e3)
Noisy1 = observed_counts[["counts"]]


##1b. Simulating Noisycounts using default parameters
set.seed(111)
data(gene_len_pool)
gene_len <- sample(gene_len_pool, ngenes, replace = FALSE)

observed_counts <- True2ObservedCounts(true_counts=true_counts_res[[1]], meta_cell=true_counts_res[[3]], 
                                       protocol="nonUMI", alpha_mean=0.1, alpha_sd=0.05, gene_len=gene_len,
                                       depth_mean=1e5, depth_sd=3e3)
Noisy1 = observed_counts[["counts"]]





 ## 2a. Simulating Noisy Counts (UMI)

set.seed(1)
data(gene_len_pool)
gene_len <- sample(gene_len_pool, ngenes, replace = FALSE)
m = 0.009
alpha_mean = rnorm(1, mean = m, sd = 0.001)
d = 85000
depht_mean = round(rnorm(1, mean = d, sd = d*0.20))

observed_counts <- True2ObservedCounts(true_counts=true_counts_res[[1]], meta_cell=true_counts_res[[3]], 
                                       protocol="UMI",  alpha_mean = alpha_mean ,  alpha_sd = alpha_mean*0.25, gene_len=gene_len,
                                       depth_mean = depht_mean, nPCR1 = 14, depth_sd = depht_mean*0.25)
Noisy1 = observed_counts[["counts"]]



## 2b. Simulating Noisy Counts (nonUMI)

set.seed(1)
data(gene_len_pool)
gene_len <- sample(gene_len_pool, ngenes, replace = FALSE)
m = 0.009
alpha_mean = rnorm(1, mean = m, sd = 0.001)
d = 85000
depht_mean = round(rnorm(1, mean = d, sd = d*0.20))

observed_counts <- True2ObservedCounts(true_counts=true_counts_res[[1]], meta_cell=true_counts_res[[3]], 
                                       protocol="nonUMI",  alpha_mean = alpha_mean ,  alpha_sd = alpha_mean*0.25, gene_len=gene_len,
                                       depth_mean = depht_mean, nPCR1 = 14, depth_sd = depht_mean*0.25)
Noisy2 = observed_counts[["counts"]]


tmp_so <- summarizedToSeurat(observed_counts)
clustering_result = calculate_clustering(tmp_so)





## 3a. Simulating Noisy Counts (UMI) -  (10000by3000)
set.seed(1)
data(gene_len_pool)
gene_len <- sample(gene_len_pool, ngenes, replace = FALSE)
m = 0.009
alpha_mean = rnorm(1, mean = m, sd = 0.001)
d = 85000
depht_mean = round(rnorm(1, mean = d, sd = d*0.20))

observed_counts <- True2ObservedCounts(true_counts=true_counts_res[[1]], meta_cell=true_counts_res[[3]], 
                                       protocol="UMI",  alpha_mean = alpha_mean ,  alpha_sd = alpha_mean*0.25, gene_len=gene_len,
                                       depth_mean = depht_mean, nPCR1 = 14, depth_sd = depht_mean*0.25)
Noisy1 = observed_counts[["counts"]]




## 3b. Simulating Noisy Counts (nonUMI)

set.seed(1)
data(gene_len_pool)
gene_len <- sample(gene_len_pool, ngenes, replace = FALSE)
m = 0.009
alpha_mean = rnorm(1, mean = m, sd = 0.001)
d = 85000
depht_mean = round(rnorm(1, mean = d, sd = d*0.20))

observed_counts <- True2ObservedCounts(true_counts=true_counts_res[[1]], meta_cell=true_counts_res[[3]], 
                                       protocol="nonUMI",  alpha_mean = alpha_mean ,  alpha_sd = alpha_mean*0.25, gene_len=gene_len,
                                       depth_mean = depht_mean, nPCR1 = 14, depth_sd = depht_mean*0.25)
Noisy2 = observed_counts[["counts"]]


# Visualizing the clusters
tmp_so <- summarizedToSeurat(observed_counts)
clustering_result = calculate_clustering(tmp_so)
  
