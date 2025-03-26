# WAVDeSc: Wavelet Denoising for single-cell RNA-seq Data

WAVDeSc is a matlab tool for denoising single-cell RNA-seq. WAVDesC leverages biorthogonal wavelet transforms to decompose scRNA-seq data into different frequency components and applies Bayesian thresholding to remove noise. The pipeline comprises three main phases: Signal Decomposition, Thresholding, and Signal Reconstruction aimed at prducing a denoised output

## Subheader
The pipeline depends on the Wavelet Toolbox.
- Check if it is installed:</br>
  `licence('test', 'Wavelet_Toolbox')`
