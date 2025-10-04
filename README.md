<p align="center">
  <strong><h1>WAVDeSc: Wavelet Denoising for single-cell RNA-seq Data</h1></strong>
</p>

This repository contains the data, scripts and notebooks used to generate the figures in our analysis of **scRNA-seq denoising methods** (DCA, ENHANCE, MAGIC, SAVER) and **WAVDESC** the tool we propose. WAVDeSc leverages biorthogonal wavelet transforms to decompose scRNA-seq data into different frequency components and applies Bayesian thresholding to remove noise. The pipeline comprises three main phases: Signal Decomposition, Thresholding, and Signal Reconstruction aimed at producing a denoised output. The approach enables the recovery of technical zeros and enhances quality while preserving important biological signals and improving downstream analyses.  

All scripts are stored in the [`scripts/`](./scripts) directory.

---

## Repository Structure
├── Denoising Methods/ </br>
├── Evaluation methods codes/ </br>
├── datasets/ </br>
├── results/ </br>
├── scripts/ </br>
└── README.md 



---

## Figures and Corresponding Scripts

| **Figure** | **Description** | **Script** |
|------------|-----------------|------------|
| Figure 3   | Heatmaps | [`scripts/notebooks/fig_3/Heatmaps.ipynb`](scripts/notebooks/fig_3/Heatmaps.ipynb) |
| Figure 4   | Evaluation Methods| [`scripts/DCA_tsne.ipynb`](./scripts/DCA_tsne.ipynb) |


---

## Usage

1. Clone the repository  
   ```bash
   git clone https://github.com/<your-username>/<repo-name>.git
   cd <repo-name>

## Features  

- Reads scRNA-seq data  
- Uses wavelet transform for noise reduction  

## Dependencies  

Ensure that you have MATLAB installed with:  
  - Wavelet Toolbox:</br>  
  `licence('test', 'Wavelet_Toolbox')`  

## License  

This project is open-source.  
