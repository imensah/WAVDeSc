### We simulated scRNA-seq data with different dimensions of rows and columns
#using the symsim software. Firstly we simulated the groundtruth to represent
#the truth data and then visualised the results using clusters.
#We then added on noise to the simulated groundtruth to produce the
#noisy data.

# Load necessary libraries
libraries <- c(
  "SymSim", "stats", "ggplot2", "cowplot", "SeuratObject", "SummarizedExperiment",
  "fdm2id", "coop", "Metrics", "ggpubr", "Seurat", "doParallel", "foreach", "reshape2"
)
# Load all libraries and suppress messages
invisible(lapply(libraries, function(lib) suppressMessages(library(lib, character.only = TRUE))))

# Set up parallel processing
cores <- detectCores()-1  # Detect number of cores available for parallel processing, and reserve one.
cl <- makeCluster(cores) 
registerDoParallel(cl) 

# Specify the simulation dimension
ngenes <<- 100
ncells <<- 50

work_dir <- "/home/jay/Desktop/WAVDeSc_repo/Data Simulations"

# Simulating True Counts based on simulation dimension
simulate_truecount <- function(ngenes, ncells) {
  phyla1 <- Phyla5() # First phylogenetic tree
  phyla2 <- read.tree(system.file("extdata", "Newick_ABCDE.txt", package = "SymSim")) # Second tree
  
  # Select simulation parameters based on dimension(ngenes and ncells)
  gt_params <- if (ngenes <= 1000 & ncells<= 500) {
    list(min_popsize = 10, i_minpop = 3, Sigma = 0.4, randseed = 10, phyla = phyla1)
  } else {
    list(min_popsize = 100, i_minpop = 30, Sigma = 0.5, randseed = 10, phyla = phyla2)
  }
  
  # Generate true counts using SimulateTrueCounts function 
  true_counts_res <<- SymSim::SimulateTrueCounts(
    ncells_total = ncells, 
    ngenes = ngenes,
    min_popsize = gt_params$min_popsize,
    i_minpop = gt_params$i_minpop,
    nevf = 8, 
    evf_type = "discrete", 
    n_de_evf = 5, 
    vary = "s", 
    Sigma = gt_params$Sigma, 
    phyla = gt_params$phyla,
    randseed = gt_params$randseed
  )
}

# Extraction and checking the true count matrix
true_counts_res <- simulate_truecount(ngenes,ncells)
groundtruth <- true_counts_res[["counts"]]
df_title <- paste("groundtruth (", ngenes, " x ", ncells, ")", sep = "") #setting title for the Dataframe
write.csv(groundtruth, file = paste0(work_dir, "/",df_title,".csv"), row.names = FALSE)



# Simulate noisy counts based on protocol (UMI and nonUMI)
simulate_noisycount <- function(protocol) {
  set.seed(111) # Set a seed for reproducibility
  data(gene_len_pool) # Load gene length pool
  gene_len <- sample(gene_len_pool, ngenes, replace = FALSE)
  
  # Set protocol-specific parameters (UMI or nonUMI)
  n_params <- if (protocol == "UMI") {
    list(alpha_mean = 0.1, alpha_sd = 0.05, depth_mean = 1e5, depth_sd = 3e3)
  } else {
    m <- 0.009  
    d <- 85000  
    list(
      alpha_mean = rnorm(1, mean = m, sd = 0.001),  
      alpha_sd = m * 0.25,  
      depth_mean = round(rnorm(1, mean = d, sd = d * 0.20)),  
      depth_sd = d * 0.25  
    )
  }
  


  # Generate noisy counts using the True2ObservedCounts function
  observed_counts <- True2ObservedCounts(
    true_counts = true_counts_res[[1]],
    meta_cell = true_counts_res[[3]],
    protocol = protocol,
    gene_len = gene_len,
    alpha_mean = n_params$alpha_mean,
    alpha_sd = n_params$alpha_sd, 
    depth_mean = n_params$depth_mean,
    depth_sd = n_params$depth_sd,
    nPCR1 = 14
  )
}


## Extraction and checking the noisy count matrix for both UMI and nonUMI
# Run simulation for noisy counts (UMI) and view results
observed_counts_umi <- simulate_noisycount("UMI")
noisy_counts_umi <- observed_counts_umi[["counts"]]
df_title_umi <- paste("noisy counts UMI (", ngenes, " x ", ncells, ")", sep = "") #setting title for the Dataframe
write.csv(noisy_counts_umi, file = paste0(work_dir, "/",df_title_umi,".csv"), row.names = FALSE)

# Run simulation for noisy counts (nonUMI) and view results
observed_counts_nonumi <- simulate_noisycount("nonUMI")
noisy_counts_nonumi <- observed_counts_nonumi[["counts"]]
df_title_nonumi <- paste("noisy counts nonUMI (", ngenes, " x ", ncells, ")", sep = "") #setting title for the Dataframe
write.csv(noisy_counts_nonumi, file = paste0(work_dir, "/",df_title_nonumi,".csv"), row.names = FALSE)


# Plotting clusters for Various Simulations
# Plot t-SNE for true counts 
tsne_true_counts <- PlotTsne(
  meta = true_counts_res[[3]],
  data = log2(true_counts_res[[1]] + 1),
  evf_type = "discrete",
  n_pc = 10,
  label = 'pop',
  saving = FALSE,
  perplexity = 2,
  plotname = paste("discrete populations (true counts)(", ngenes, " x ", ncells, ")", sep = "")
)

# Plot t-SNE for noisy counts (UMI)
tsne_noisy_counts_umi <- PlotTsne(
  meta = true_counts_res[[3]],
  data = log2(noisy_counts_umi + 1),  
  evf_type = "discrete",
  n_pc = 10,
  label = 'pop',
  saving = FALSE,
  perplexity = 2,
  plotname = paste("Noisy Data UMI: (", ngenes, " x ", ncells, ")", sep = "")
)

# Plot t-SNE for noisy counts (non-UMI)
tsne_noisy_counts_nonumi <- PlotTsne(
  meta = true_counts_res[[3]],
  data = log2(noisy_counts_nonumi + 1),  
  evf_type = "discrete",
  n_pc = 10,
  label = 'pop',
  saving = FALSE,
  perplexity = 2,
  plotname = paste("Noisy Data nonUMI: (", ngenes, " x ", ncells, ")", sep = "")
)

# Combine all t-SNE plots into one figure
combined_tsne_plot <- cowplot::plot_grid(
  tsne_true_counts[[2]], 
  tsne_noisy_counts_umi[[2]],
  tsne_noisy_counts_nonumi[[2]],
  ncol = 3,  # Arrange plots in 3 columns
  align = "h"  # Align plots horizontally           
)

# Print the combined t-SNE plot
print(combined_tsne_plot)


# Stop the cluster
stopCluster(cl)