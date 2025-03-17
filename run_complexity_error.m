%% this script runs the main experment to analyse the algorithmic complexity of dft vs fft and relative error

%% Add functions folder 'src' to the MATLAB search path
addpath(fullfile('src'));

%% ===SETUP AND PARAMETERS===
% Sampling frequency (Hz)
real_rate = 300;
% Input frequency for the ideal sine wave (Hz) - modify as needed
input_freq = 5;

% Define logarithmic range of sequence lengths (using powers of 2)
k_min = 5;    % smallest exponent (e.g., N = 2^5 = 32)
k_max = 14;   % largest exponent (e.g., N = 2^14 = 16384)
N_values = 2.^(k_min:k_max);  % sequence lengths

% Number of times to run each experiment (enter desired number of runs)
% This is so we can compute error bars
num_expt = 5;  % Change this number as needed

% Preallocate matrices to store runtime measurements and RMSE values.
% Rows correspond to different sequence lengths; columns correspond to experiment repetitions.
runtime_dft_all = zeros(length(N_values), num_expt);
runtime_fft_all = zeros(length(N_values), num_expt);
rmse_dft_all    = zeros(length(N_values), num_expt);
rmse_fft_all    = zeros(length(N_values), num_expt);

%% ===ALGORITHMIC COMPLEXITY AND ERROR ANALYSIS===
% Loop over each sequence length
for idx = 1:length(N_values)
    N = N_values(idx);
    
    % Create a time vector for N samples
    t = (0:N-1)' / real_rate;
    
    % Generate the ideal sine wave signal
    signal = sin(2*pi*input_freq*t);
    
    % Run each experiment multiple times for the current sequence length
    for expt = 1:num_expt
        % Measure runtime for the custom dft_vectorized function
        tic;
        dft_signal = dft_vectorized(signal);
        runtime_dft_all(idx, expt) = toc;
        
        % Measure runtime for the custom fft_vectorized function
        tic;
        fft_signal = fft_vectorized(signal);
        runtime_fft_all(idx, expt) = toc;

        % --- Compute reconstruction using MATLAB's IDFT (ifft) ---
        rec_dft = ifft(dft_signal);
        rec_fft = ifft(fft_signal);

        % --- Compute RMSE between original signal and reconstructed signal ---
        % Since ifft gives us complex values, we need to take the absolute value of the difference first before doing rmse
        % the . dot operator tells matlab to perform element wise operations
        rmse_dft_all(idx, expt) = sqrt(mean(abs(signal - rec_dft).^2));
        rmse_fft_all(idx, expt) = sqrt(mean(abs(signal - rec_fft).^2));
    end
end

% Compute mean and standard deviation for each sequence length (across runs)
mean_runtime_dft = mean(runtime_dft_all, 2); %calculate mean across rows
std_runtime_dft  = std(runtime_dft_all, 0, 2); %std across rows
mean_runtime_fft = mean(runtime_fft_all, 2);
std_runtime_fft  = std(runtime_fft_all, 0, 2);

% Compute mean and standard deviation of RMSE across experiments for each sequence length
mean_rmse_dft = mean(rmse_dft_all, 2);
std_rmse_dft  = std(rmse_dft_all, 0, 2);
mean_rmse_fft = mean(rmse_fft_all, 2);
std_rmse_fft  = std(rmse_fft_all, 0, 2);

%% PLOTTING THE COMPLEXITY RESULTS WITH ERROR BARS
figure;
% Plot DFT mean runtime with error bars
errorbar(N_values, mean_runtime_dft, std_runtime_dft, '-o', 'LineWidth', 2);
hold on;
% Plot FFT mean runtime with error bars
errorbar(N_values, mean_runtime_fft, std_runtime_fft, '-s', 'LineWidth', 2);
xlabel('Sequence Length (N)');
ylabel('Runtime (seconds)');
title('Algorithmic Complexity: dft\_vectorized vs fft\_vectorized');
legend('dft\_vectorized', 'fft\_vectorized', 'Location', 'NorthWest');
grid on;

% Set axes to logarithmic scale for both X and Y axes
set(gca, 'XScale', 'log', 'YScale', 'log');
hold off;

% Save the plot as a PNG image with specified renderer settings
print('-dpng', '-r300', '-painters', 'complexity.png');
fprintf('Algorithm complexity plot saved\n');

%% PLOTTING RMSE (RECONSTRUCTION ERROR) WITH ERROR BARS
figure;
% Plot RMSE for reconstruction using dft_vectorized with MATLAB's ifft
errorbar(N_values, mean_rmse_dft, std_rmse_dft, '-o', 'LineWidth', 2);
hold on;
% Plot RMSE for reconstruction using fft_vectorized with MATLAB's ifft
errorbar(N_values, mean_rmse_fft, std_rmse_fft, '-s', 'LineWidth', 2);
xlabel('Sequence Length (N)');
ylabel('RMSE');
title('Reconstruction Error (RMSE) via MATLAB IDFT');
legend('DFT-based reconstruction', 'FFT-based reconstruction', 'Location', 'NorthWest');
grid on;
set(gca, 'XScale', 'log'); % Logarithmic x-axis is usually sufficient here
hold off;

% Save the RMSE plot as a PNG image
print('-dpng', '-r300', '-painters', 'rmse_analysis.png');
fprintf('RMSE analysis plot saved as rmse_analysis.png\n');
