### We simulated scRNA-seq data with different dimensions of rows and columns
#using the symsim software. Firstly we simulated the groundtruth to represent
#the truth data and then visualised the results using clusters.
#We then used then added on noise to the simulated groundtruth to produce the
#noisy data.

## Load all packages needed for the analysis
# List of libraries to load
libraries <- c(
  "SymSim",
  "stats",
  "coop",
  "ggplot2",
  "cowplot", 
  "fdm2id", 
  "Metrics", 
  "ggpubr", 
  "Seurat", 
  "SeuratObject", 
  "SummarizedExperiment")

# Load all libraries using lapply
lapply(libraries, library, character.only = TRUE)



## NB: Here simulations are done with some of the the default simulation
#parameters used in the symsim paper.

simulate_truecount <- function(ncells,ngenes) {
    phyla1 <- Phyla5()
    phyla2 <- read.tree(system.file("extdata", "Newick_AB.txt", package = "SymSim"))
    phyla3 <- read.tree(system.file("extdata", "Newick_ABCDE.txt", package = "SymSim"))

## 1. Simulating Truecount with dimensions 500 genes by 300 cells with
# five(5) cell populations
    if (ngenes <= 500 & ncells <= 300) {
        min_popsize <- 50
        i_minpop <- 2
        nevf <- 10
        Sigma <- 0.4
        phyla = phyla1
        randseed = 0
    }
## 2. Simulating Truecount with dimensions 6000 genes by 3000 cells
# with three(3) cell populations
    else if (ngenes <= 6000 & ncells <= 3000) {
        plot(phyla2)

        min_popsize <- 800
        i_minpop <- 3
        nevf <- 20
        gene_effect_prob <- 0.1
        gene_effects_sd <- 0.5
        Sigma <- 0.5
        phyla <- phyla2
        randseed <- 10
    }
## 3. Simulating Truecount with dimensions 10000 genes by 3000 cells
# with five cell populations
    else {
        min_popsize <- 600
        nevf <- 60
        gene_effect_prob <- 0.1
        gene_effects_sd <- 0.5
        Sigma <- 0.5
        phyla <- phyla1
        randseed <- 10
    }
    true_counts_res <- SimulateTrueCounts(
        ncells_total = ncells,
        min_popsize = min_popsize,
        i_minpop = ifelse(exists("i_minpop"),i_minpop,NULL),
        ngenes = ngenes,
        nevf = nevf,
        evf_type = "discrete",
        n_de_evf = 9,
        vary = "s",
        gene_effect_prob = ifelse(exists("gene_effect_prob"), gene_effect_prob,NULL),
        gene_effects_sd = ifelse(exists("gene_effect_sd"), gene_effect_prob, NULL),
        Sigma = Sigma,
        phyla = phyla,
        randseed = randseed
    )

#Extraction and checking the true count matrix for ground truth two
    groundtruth <- true_counts_res[["counts"]]
    View(groundtruth)
}

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