function [denoised_data, metrics, params_used] = WAVDeSc(input_data, varargin)
% WAVDeSc - Wavelet-Based Denoising for Single-Cell RNA Sequencing Data
%
% DESCRIPTION:
%   WAVDeSc applies wavelet-based denoising to scRNA-seq data using discrete
%   wavelet transform (DWT) with Bayesian thresholding to reduce technical
%   noise while preserving biological signal.
%
% SYNTAX:
%   denoised_data = WAVDeSc(input_data)
%   denoised_data = WAVDeSc(input_data, 'Parameter', Value, ...)
%   [denoised_data, metrics, params_used] = WAVDeSc(...)
%
% INPUTS:
%   input_data - Input can be:
%                1) Numeric matrix (genes × cells or cells × genes)
%                2) Table with gene names in first column/row
%                3) CSV file path (string or char)
%
% OPTIONAL PARAMETERS (Name-Value pairs):
%   'Wavelet'          - Wavelet function (default: 'db6')
%                        Options: 'db4', 'db6', 'db8', 'bior2.6', etc.
%   'Orientation'      - Data orientation (default: 'auto')
%                        'auto': auto-detect, 'genes_rows': genes as rows,
%                        'genes_cols': genes as columns
%   'DenoisingMethod'  - Denoising method (default: 'Bayes')
%                        Options: 'Bayes', 'UniversalThreshold', 'Minimax'
%   'ThresholdRule'    - Thresholding rule (default: 'Hard')
%                        Options: 'Hard', 'Soft'
%   'NoiseEstimate'    - Noise estimation (default: 'LevelDependent')
%                        Options: 'LevelDependent', 'LevelIndependent'
%   'DecompositionLevel' - Decomposition level (default: 'auto')
%                        'auto': automatically calculated, or specify integer
%   'SaveOutput'       - Save denoised data (default: false)
%   'OutputPath'       - Path to save output (default: './WAVDeSc_output.csv')
%   'ComputeMetrics'   - Compute performance metrics if ground truth provided
%                        (default: false)
%   'GroundTruth'      - Ground truth data for metric calculation
%   'Verbose'          - Display progress messages (default: true)
%   'PlotResults'      - Generate visualization plots (default: false)
%
% OUTPUTS:
%   denoised_data - Denoised scRNA-seq matrix (same orientation as input)
%   metrics       - Structure containing performance metrics (if computed)
%   params_used   - Structure containing all parameters used
%
% EXAMPLES:
%   % Example 1: Basic usage with file path
%   denoised = WAVDeSc('path/to/data.csv');
%
%   % Example 2: With custom wavelet and parameters
%   denoised = WAVDeSc(data_matrix, 'Wavelet', 'db8', 'ThresholdRule', 'Soft');
%
%   % Example 3: With ground truth for evaluation
%   [denoised, metrics] = WAVDeSc(noisy_data, 'GroundTruth', true_data, ...
%                                  'ComputeMetrics', true, 'PlotResults', true);
%
%   % Example 4: Save output with custom path
%   denoised = WAVDeSc('data.csv', 'SaveOutput', true, ...
%                      'OutputPath', 'results/denoised_data.csv');
%
% REFERENCES:
%   [1] Your paper citation here
%   [2] Wavelet toolbox documentation
%
% AUTHORS:
%   Isabel Mensah, Justice Kwame Appati, Samson Pandam Salifu, 
%   Peter Amoako-Yirenkyi
%
% VERSION: 1.0
% DATE: 2025

%% Input parsing and validation
p = inputParser;
addRequired(p, 'input_data');
addParameter(p, 'Wavelet', 'db6', @(x) ischar(x) || isstring(x));
addParameter(p, 'Orientation', 'auto', @(x) any(validatestring(x, {'auto', 'genes_rows', 'genes_cols'})));
addParameter(p, 'DenoisingMethod', 'Bayes', @(x) any(validatestring(x, {'Bayes', 'UniversalThreshold', 'Minimax'})));
addParameter(p, 'ThresholdRule', 'Hard', @(x) any(validatestring(x, {'Hard', 'Soft'})));
addParameter(p, 'NoiseEstimate', 'LevelDependent', @(x) any(validatestring(x, {'LevelDependent', 'LevelIndependent'})));
addParameter(p, 'DecompositionLevel', 'auto');
addParameter(p, 'SaveOutput', false, @islogical);
addParameter(p, 'OutputPath', './WAVDeSc_output.csv', @(x) ischar(x) || isstring(x));
addParameter(p, 'ComputeMetrics', false, @islogical);
addParameter(p, 'GroundTruth', []);
addParameter(p, 'Verbose', true, @islogical);
addParameter(p, 'PlotResults', false, @islogical);

parse(p, input_data, varargin{:});
params = p.Results;

%% Load and prepare data
if params.Verbose
    fprintf('\n========================================\n');
    fprintf('WAVDeSc: Wavelet-Based scRNA-seq Denoising\n');
    fprintf('========================================\n\n');
    fprintf('Step 1/5: Loading and preparing data...\n');
end

[data_matrix, gene_names, cell_names, original_orientation] = load_scrna_data(params.input_data, params.Verbose);

% Store original dimensions
[n_rows, n_cols] = size(data_matrix);

% Determine data orientation
if strcmp(params.Orientation, 'auto')
    % Assume genes > cells typically
    if n_rows > n_cols
        orientation = 'genes_rows';
    else
        orientation = 'genes_cols';
    end
    if params.Verbose
        fprintf('   Auto-detected orientation: %s\n', orientation);
    end
else
    orientation = params.Orientation;
end

% Ensure data is in cells × genes format for processing
if strcmp(orientation, 'genes_rows')
    data_matrix = data_matrix';
    processing_transposed = true;
else
    processing_transposed = false;
end

[n_cells, n_genes] = size(data_matrix);

if params.Verbose
    fprintf('   Data dimensions: %d cells × %d genes\n', n_cells, n_genes);
    fprintf('   Sparsity: %.2f%%\n', 100 * sum(data_matrix(:) == 0) / numel(data_matrix));
end

%% Determine decomposition level
if params.Verbose
    fprintf('\nStep 2/5: Calculating decomposition level...\n');
end

if strcmp(params.DecompositionLevel, 'auto')
    declev = [n_cells, n_genes];
    declev_d = wmaxlev(declev, params.Wavelet);
    if params.Verbose
        fprintf('   Auto-calculated decomposition level: %d\n', declev_d);
    end
else
    declev_d = params.DecompositionLevel;
    if params.Verbose
        fprintf('   User-specified decomposition level: %d\n', declev_d);
    end
end

%% Perform wavelet denoising
if params.Verbose
    fprintf('\nStep 3/5: Performing wavelet denoising...\n');
    fprintf('   Wavelet: %s\n', params.Wavelet);
    fprintf('   Denoising method: %s\n', params.DenoisingMethod);
    fprintf('   Threshold rule: %s\n', params.ThresholdRule);
    fprintf('   Noise estimate: %s\n', params.NoiseEstimate);
    tic;
end

try
    denoised_matrix = wdenoise(data_matrix, declev_d, ...
        'Wavelet', params.Wavelet, ...
        'DenoisingMethod', params.DenoisingMethod, ...
        'ThresholdRule', params.ThresholdRule, ...
        'NoiseEstimate', params.NoiseEstimate);
catch ME
    error('WAVDeSc:DenoisingError', ...
        'Wavelet denoising failed: %s\nTry reducing decomposition level or changing wavelet.', ...
        ME.message);
end

if params.Verbose
    elapsed = toc;
    fprintf('   Denoising completed in %.2f seconds\n', elapsed);
end

%% Post-processing
if params.Verbose
    fprintf('\nStep 4/5: Post-processing...\n');
end

% Replace negative values with zero
n_negative = sum(denoised_matrix(:) < 0);
denoised_matrix(denoised_matrix < 0) = 0;

if params.Verbose
    fprintf('   Negative values replaced with zero: %d (%.2f%%)\n', ...
        n_negative, 100 * n_negative / numel(denoised_matrix));
    fprintf('   Final sparsity: %.2f%%\n', ...
        100 * sum(denoised_matrix(:) == 0) / numel(denoised_matrix));
end

% Transpose back to original orientation if needed
if processing_transposed
    denoised_data = denoised_matrix';
else
    denoised_data = denoised_matrix;
end

%% Compute metrics if requested
metrics = struct();
if params.ComputeMetrics && ~isempty(params.GroundTruth)
    if params.Verbose
        fprintf('\nStep 5/5: Computing performance metrics...\n');
    end
    
    % Ensure ground truth has same orientation
    ground_truth = params.GroundTruth;
    if processing_transposed
        if size(ground_truth, 1) == n_genes && size(ground_truth, 2) == n_cells
            ground_truth = ground_truth';
        end
    end
    
    % Compute metrics
    metrics = compute_metrics(denoised_matrix, ground_truth, data_matrix, params.Verbose);
    
elseif params.Verbose
    fprintf('\nStep 5/5: Skipping metrics computation (no ground truth provided)\n');
end

%% Save output
if params.SaveOutput
    if params.Verbose
        fprintf('\nSaving denoised data to: %s\n', params.OutputPath);
    end
    
    % Create table with gene and cell names if available
    if ~isempty(gene_names) && ~isempty(cell_names)
        if strcmp(original_orientation, 'genes_rows')
            output_table = array2table(denoised_data, ...
                'RowNames', gene_names, 'VariableNames', cell_names);
        else
            output_table = array2table(denoised_data, ...
                'RowNames', cell_names, 'VariableNames', gene_names);
        end
        writetable(output_table, params.OutputPath, 'WriteRowNames', true);
    else
        dlmwrite(params.OutputPath, denoised_data, 'precision', '%.6f');
    end
    
    if params.Verbose
        fprintf('   Output saved successfully!\n');
    end
end

%% Generate plots
if params.PlotResults
    generate_plots(data_matrix, denoised_matrix, metrics, params);
end

%% Prepare output parameters
params_used = struct();
params_used.Wavelet = params.Wavelet;
params_used.DecompositionLevel = declev_d;
params_used.DenoisingMethod = params.DenoisingMethod;
params_used.ThresholdRule = params.ThresholdRule;
params_used.NoiseEstimate = params.NoiseEstimate;
params_used.InputDimensions = [n_rows, n_cols];
params_used.ProcessingDimensions = [n_cells, n_genes];
params_used.Orientation = orientation;

if params.Verbose
    fprintf('\n========================================\n');
    fprintf('WAVDeSc denoising completed successfully!\n');
    fprintf('========================================\n\n');
end

end

%% Helper Functions

function [data_matrix, gene_names, cell_names, orientation] = load_scrna_data(input_data, verbose)
% Load scRNA-seq data from various input formats

gene_names = {};
cell_names = {};
orientation = 'unknown';

if ischar(input_data) || isstring(input_data)
    % Input is a file path
    if verbose
        fprintf('   Loading data from file: %s\n', input_data);
    end
    
    try
        data_table = readtable(input_data);
        
        % Check if first column/row contains names
        if iscell(data_table{:,1}) || isstring(data_table{:,1})
            gene_names = data_table{:,1};
            data_table(:,1) = [];
        end
        
        % Get cell names from variable names
        cell_names = data_table.Properties.VariableNames;
        
        % Convert to numeric matrix
        data_matrix = table2array(data_table);
        orientation = 'genes_rows';
        
    catch
        % Try as plain numeric file
        data_matrix = dlmread(input_data);
    end
    
elseif istable(input_data)
    % Input is a table
    if iscell(input_data{:,1}) || isstring(input_data{:,1})
        gene_names = input_data{:,1};
        input_data(:,1) = [];
    end
    cell_names = input_data.Properties.VariableNames;
    data_matrix = table2array(input_data);
    orientation = 'genes_rows';
    
elseif isnumeric(input_data)
    % Input is already a numeric matrix
    data_matrix = input_data;
    
else
    error('WAVDeSc:InvalidInput', ...
        'Input must be a numeric matrix, table, or file path (string/char)');
end

% Validate data
if ~isnumeric(data_matrix)
    error('WAVDeSc:InvalidData', 'Data must be numeric');
end

if any(isnan(data_matrix(:)))
    warning('WAVDeSc:NaNValues', ...
        'Data contains NaN values. Replacing with zeros.');
    data_matrix(isnan(data_matrix)) = 0;
end

if any(isinf(data_matrix(:)))
    error('WAVDeSc:InfValues', 'Data contains Inf values');
end

end

function metrics = compute_metrics(denoised, ground_truth, noisy, verbose)
% Compute performance metrics comparing denoised data to ground truth

metrics = struct();

% Ensure same dimensions
if ~isequal(size(denoised), size(ground_truth))
    error('WAVDeSc:DimensionMismatch', ...
        'Denoised data and ground truth must have same dimensions');
end

% RMSE
metrics.RMSE = sqrt(mean((ground_truth(:) - denoised(:)).^2));
metrics.RMSE_noisy = sqrt(mean((ground_truth(:) - noisy(:)).^2));

% Correlation
metrics.Correlation = corr(ground_truth(:), denoised(:), 'type', 'Pearson');
metrics.Correlation_Spearman = corr(ground_truth(:), denoised(:), 'type', 'Spearman');

% MSE
metrics.MSE = mean((ground_truth(:) - denoised(:)).^2);

% MAE
metrics.MAE = mean(abs(ground_truth(:) - denoised(:)));

% R-squared
SS_res = sum((ground_truth(:) - denoised(:)).^2);
SS_tot = sum((ground_truth(:) - mean(ground_truth(:))).^2);
metrics.R_squared = 1 - (SS_res / SS_tot);

% Signal-to-Noise Ratio improvement
signal_power = var(ground_truth(:));
noise_power_before = var(ground_truth(:) - noisy(:));
noise_power_after = var(ground_truth(:) - denoised(:));
metrics.SNR_improvement_dB = 10 * log10(noise_power_before / noise_power_after);

% Sparsity metrics
metrics.Sparsity_truth = 100 * sum(ground_truth(:) == 0) / numel(ground_truth);
metrics.Sparsity_noisy = 100 * sum(noisy(:) == 0) / numel(noisy);
metrics.Sparsity_denoised = 100 * sum(denoised(:) == 0) / numel(denoised);

if verbose
    fprintf('   RMSE: %.4f (noisy: %.4f)\n', metrics.RMSE, metrics.RMSE_noisy);
    fprintf('   Pearson correlation: %.4f\n', metrics.Correlation);
    fprintf('   R-squared: %.4f\n', metrics.R_squared);
    fprintf('   SNR improvement: %.2f dB\n', metrics.SNR_improvement_dB);
end

end

function generate_plots(noisy, denoised, metrics, params)
% Generate visualization plots

figure('Name', 'WAVDeSc Results', 'Position', [100, 100, 1200, 800]);

% Plot 1: Sample gene expression (first gene)
subplot(2, 3, 1);
plot(noisy(:, 1), 'b-', 'LineWidth', 1.5);
hold on;
plot(denoised(:, 1), 'r-', 'LineWidth', 1.5);
if ~isempty(params.GroundTruth)
    if size(params.GroundTruth, 2) == size(noisy, 2)
        plot(params.GroundTruth(:, 1), 'g--', 'LineWidth', 1.5);
        legend('Noisy', 'Denoised', 'Ground Truth');
    else
        legend('Noisy', 'Denoised');
    end
else
    legend('Noisy', 'Denoised');
end
title('Sample Gene Expression (Gene 1)');
xlabel('Cells');
ylabel('Expression');
grid on;

% Plot 2: Sparsity comparison
subplot(2, 3, 2);
sparsity_noisy = 100 * sum(noisy == 0, 1) / size(noisy, 1);
sparsity_denoised = 100 * sum(denoised == 0, 1) / size(denoised, 1);
plot(sparsity_noisy, 'b-', 'LineWidth', 1.5);
hold on;
plot(sparsity_denoised, 'r-', 'LineWidth', 1.5);
legend('Noisy', 'Denoised');
title('Sparsity per Gene');
xlabel('Genes');
ylabel('Sparsity (%)');
grid on;

% Plot 3: Distribution comparison
subplot(2, 3, 3);
histogram(log1p(noisy(:)), 50, 'Normalization', 'probability', 'FaceAlpha', 0.5);
hold on;
histogram(log1p(denoised(:)), 50, 'Normalization', 'probability', 'FaceAlpha', 0.5);
legend('Noisy', 'Denoised');
title('Expression Distribution (log1p)');
xlabel('log1p(Expression)');
ylabel('Probability');
grid on;

% Plot 4: Heatmap comparison (sample)
n_sample = min(50, size(noisy, 1));
genes_sample = min(100, size(noisy, 2));

subplot(2, 3, 4);
imagesc(noisy(1:n_sample, 1:genes_sample));
colorbar;
title('Noisy Data (Sample)');
xlabel('Genes');
ylabel('Cells');

subplot(2, 3, 5);
imagesc(denoised(1:n_sample, 1:genes_sample));
colorbar;
title('Denoised Data (Sample)');
xlabel('Genes');
ylabel('Cells');

% Plot 6: Metrics (if available)
subplot(2, 3, 6);
if ~isempty(fieldnames(metrics))
    metric_names = {'RMSE', 'Correlation', 'R_squared'};
    metric_values = [metrics.RMSE, metrics.Correlation, metrics.R_squared];
    bar(metric_values);
    set(gca, 'XTickLabel', metric_names);
    title('Performance Metrics');
    ylabel('Value');
    grid on;
    
    % Add text annotations
    for i = 1:length(metric_values)
        text(i, metric_values(i), sprintf('%.3f', metric_values(i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
else
    text(0.5, 0.5, 'No metrics available', ...
        'HorizontalAlignment', 'center', 'FontSize', 12);
    axis off;
end

sgtitle(sprintf('WAVDeSc Results - Wavelet: %s, Method: %s', ...
    params.Wavelet, params.DenoisingMethod));

end
