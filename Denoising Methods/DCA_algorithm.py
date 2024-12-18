import pandas as pd
import os
import re
import numpy as np
import scanpy as sc
import dca.api as d
from anndata import AnnData

# Set the path to the folder containing the noisy datasets
data_directory = "/home/jay/Desktop/WAVDeSc_repo/Data Simulations"

# Algorithm for denoising
algorithm = "DCA"
inputDir = f"{data_directory}/Noisy_files/noisy counts nonUMI (100 x 50).csv"
outputDir = f"{data_directory}/denoised_{algorithm}"

# Ensure the output directory exists
if not os.path.isdir(outputDir):
    os.mkdir(outputDir)

# Get the list of all files in the input directory
noisy_files = os.listdir(inputDir)
# Filter the files that contain "noisy" in their name
noisy_files = [x for x in noisy_files if "noisy" in x]

# Loop through each file for denoising
for file in noisy_files:
    experiment_name = re.sub("noisy", "", file)  # Remove 'noisy' from the filename
    denoised_name = f"{outputDir}/{experiment_name}_denoised_{algorithm}.csv"  # Define the output filename for denoised data

    print(f"Processing: {denoised_name}")

    # Read the noisy dataset into an AnnData object
    adata = sc.read(f"{inputDir}/{file}")
    adata = adata.transpose()  # Transpose the data: genes in rows, cells in columns
    
    # Choose the appropriate loss function based on the dataset type
    loss = "zinb-conddisp" if "nonUMI" in experiment_name else "nb-conddisp"

    # Apply DCA denoising to the data
    adata.X = np.ceil(adata.X)  # Round the expression values (important for DCA)
    adata_denoised = d.dca(adata, copy=True, log1p=False, return_info=True, verbose=True, reduce_lr=5, scale=False, ae_type=loss)

    # Round the denoised data for consistency
    adata_denoised.X = np.round(adata_denoised.X)

    # Save the denoised data to a CSV file
    print("Saving denoised data...")
    adata_denoised_df = pd.DataFrame(adata_denoised.X, index=adata_denoised.obs_names, columns=adata_denoised.var_names)
    adata_denoised_df.to_csv(denoised_name)

    # Optionally, save size factors (if available in the data)
    if 'size_factors' in adata_denoised.obs:
        print("Saving size factors...")
        size_factors = adata_denoised.obs['size_factors']
        size_factors.to_csv(f"{outputDir}/size/{experiment_name}_dca_size.csv")

    # If it's a non-UMI dataset, save dropout probabilities
    if "nonUMI" in file:
        print("Saving dropout probabilities...")
        # Dropout probabilities are stored in 'X_dca_dropout' in the 'obsm' attribute
        dropout_data = AnnData(adata_denoised.obsm["X_dca_dropout"])
        dropout_data.obs.index = adata_denoised.obs.index
        dropout_data.var.index = adata_denoised.var.index

        # Ensure the dropout output directory exists
        os.makedirs(f"{outputDir}/dropout/", exist_ok=True)

        # Save the dropout probabilities to a CSV file
        dropout_data_df = pd.DataFrame(dropout_data.X, index=dropout_data.obs_names, columns=dropout_data.var_names)
        dropout_data_df.to_csv(f"{outputDir}/dropout/{experiment_name}_dca_dropout.csv")

    print(f"Completed processing for {experiment_name}")
