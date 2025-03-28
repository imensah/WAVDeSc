<p align="center">
<strong><span style="font-size: 500px;">WAVDeSc: Wavelet Denoising for single-cell RNA-seq Data</span></strong>
</p>

This repository contains a tool for denoising single-cell RNA-seq (scRNA-seq) data. WAVDesC leverages biorthogonal wavelet transforms to decompose scRNA-seq data into different frequency components and applies Bayesian thresholding to remove noise. The pipeline comprises three main phases: Signal Decomposition, Thresholding, and Signal Reconstruction aimed at prducing a denoised output. The approach enables the recovery of technical zeros and enhance quality while preserving important biological signals and improving downstream analyses. 

## Features

-  Reads scRNA-seq data
-  Uses wavelet transform for noise reduction
  
## Dependencies

Ensure that you have MATLAB installed with:
  - Wavelet Toolbox:</br>
  `licence('test', 'Wavelet_Toolbox')`

## License

This project is open-source.
