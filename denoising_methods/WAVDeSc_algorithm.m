clear
clc
close all


%load scRNA data
scRNA_data = readtable("path/to/scRNA-seq/transposed.csv");

% Save the first column containing cell names
cell_names = scRNA_data(:, 1);  

% Remove the first column containing cell names from the scRNA_data table
scRNA_data(:,1) = [];


%Convert the table into numeric
scRNA_data=table2array(scRNA_data);

%Estimate the decompositon level using the wavelet function
declev = [num_cells  num_genes];
declev_wf = "wavelet function";			%db4 -- appropriate wavelet function
declev_d = wmaxlev(declev,declev_wf);


%wavelet decomposition and reconstruction on noisy data using wavdec2
WAVDESC = wdenoise(scRNA_data, declev_d ,"Wavelet","wavelet function","DenoisingMethod","Bayes","ThresholdRule","Hard","NoiseEstimate","LevelDependent");


% Transpose data again
WAVDESCTT = WAVDESC';

% Replace negative expressions with zero
WAVDESCTT(WAVDESCTT < 0) = 0;

% Save results for further analysis
dlmwrite("/home/elisa/Desktop/WAVDESC_PBMC/wavdesc_pbmcD.csv", WAVDESCTT)
dlmwrite("/home/elisa/Desktop/WAVDESC_PBMC/wavdesc_pbmc.csv", WAVDESCTT)
