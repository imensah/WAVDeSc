import magic
import os
import re
import numpy as np
import scanpy as sc
import dca.api as d
from anndata import AnnData
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import resource
import time
from memory_profiler import memory_usage

# Function to get the current timestamp
def timestamp():
    return time.time()

# Function to get resource usage
def resource_usage(who):
    return resource.getrusage(who)


def find_pca_comp(adata, figName, figTitle):
    adata_scaled = adata.copy()
    sc.pp.scale(adata_scaled)
    
    # Drop rows with NaN values
    adata_scaled.X = adata_scaled.X[~np.isnan(adata_scaled.X).any(axis=1)]
    
    # Apply PCA
    pca = PCA().fit(adata_scaled.X)
    
    # Plotting the Cumulative Summation of the Explained Variance
    plt.figure()
    n = next(i for i, x in enumerate(np.cumsum(pca.explained_variance_ratio_)) if x >= 0.7)
    plt.plot(np.cumsum(pca.explained_variance_ratio_))
    plt.vlines(n, 0, 1, linestyles="dashed")
    plt.xlabel('Number of Components')
    plt.ylabel('Variance (%)')  # for each component
    plt.title(figTitle)
    plt.savefig(figName, dpi=200)
    plt.close()
    return n



# Set the path to the folder containing the noisy datasets
data_directory = "/home/jay/Desktop/WAVDeSc_repo/Data Simulations/Noisy_files"
algorithm = "MAGIC"
inputDir = f"{data_directory}/"
outputDir = f"{data_directory}/denoised_{algorithm}"


noisy_files = os.listdir(inputDir)
noisy_files = [x for x in noisy_files if "noisy counts" in x]



for files in noisy_files:
    experiment_name = re.sub("noisy counts","", files)
    denoised_name = f"{outputDir}/{experiment_name}_denoised_{algorithm}.csv"
    print(denoised_name)
    if os.path.isfile(denoised_name):
        continue

    adata = sc.read(f"{inputDir}/{files}")
    adata = adata.transpose()
    adata.X = np.expm1(adata.X)
    sc.pp.sqrt(adata)

    n = find_pca_comp(adata, 
                      figName = f"{outputDir}/figures/{experiment_name}_variance.png",
                      figTitle = f'{experiment_name} Explained Variance')

    magic_op = magic.MAGIC(t=6, n_pca = n)

    start_time, start_resources = timestamp(), resource_usage(RUSAGE_SELF)
    (mem_registered, adata_denoised) = memory_usage((magic_op.fit_transform, 
                                                    (adata,), 
                                                    {'genes': 'all_genes'}), 
                                                    retval = True, 
                                                    max_usage = True, 
                                                    include_children =True)

    end_resources, end_time = resource_usage(resource.RUSAGE_SELF), timestamp()
    real = end_time - start_time
    systime =  end_resources.ru_stime - start_resources.ru_stime
    usertime = end_resources.ru_utime - start_resources.ru_utime
    cpu_time = systime + usertime


    adata_denoised.X = adata_denoised.X**2
    print("Saving denoised_data")
    #anndata_to_csv(adata_denoised.transpose(),denoised_name)
    adata_denoised.to_df().transpose().to_csv(denoised_name)  # Save to CSV
    
    # file = open(f"{outputDir}/{algorithm}_runtime.csv", "a+")
    # file.write(f"\n{experiment_name},{algorithm},{str(cpu_time)},{str(real)},{str(mem_registered)}")
    # file.close()

    # Log runtime and memory usage
    with open(f"{outputDir}/{algorithm}_runtime.csv", "a+") as file:
        file.write(f"\n{experiment_name},{algorithm},{str(cpu_time)},{str(real)},{str(mem_registered)}")