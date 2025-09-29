<p align="center">
  <strong><h1>WAVDeSc: Wavelet Denoising for single-cell RNA-seq Data</h1></strong>
</p>

This repository contains a tool for denoising single-cell RNA-seq (scRNA-seq) data. WAVDeSc leverages biorthogonal wavelet transforms to decompose scRNA-seq data into different frequency components and applies Bayesian thresholding to remove noise. The pipeline comprises three main phases: Signal Decomposition, Thresholding, and Signal Reconstruction aimed at producing a denoised output. The approach enables the recovery of technical zeros and enhances quality while preserving important biological signals and improving downstream analyses.  

## Repo structure

## Features  

- Reads scRNA-seq data  
- Uses wavelet transform for noise reduction  

## Dependencies  

Ensure that you have MATLAB installed with:  
  - Wavelet Toolbox:</br>  
  `licence('test', 'Wavelet_Toolbox')`  

## License  

This project is open-source.  
