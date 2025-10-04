<p align="center">
  <strong><h1>WAVDeSc: Wavelet Denoising for single-cell RNA-seq Data</h1></strong>
</p>

This repository contains the data, scripts and notebooks used to generate the figures in our analysis of **scRNA-seq denoising methods** (DCA, ENHANCE, MAGIC, SAVER) and **WAVDESC** the tool we propose. WAVDeSc leverages biorthogonal wavelet transforms to decompose scRNA-seq data into different frequency components and applies Bayesian thresholding to remove noise. The pipeline comprises three main phases: Signal Decomposition, Thresholding, and Signal Reconstruction aimed at producing a denoised output. The approach enables the recovery of technical zeros and enhances quality while preserving important biological signals and improving downstream analyses.  

All scripts are stored in the [`scripts/`](./scripts) directory.

---

## Repository Structure
├── scripts/ </br>
├── data/ </br>
├── results/ </br>
└── README.md 


## Features  

- Reads scRNA-seq data  
- Uses wavelet transform for noise reduction  

## Dependencies  

Ensure that you have MATLAB installed with:  
  - Wavelet Toolbox:</br>  
  `licence('test', 'Wavelet_Toolbox')`  

## License  

This project is open-source.  
