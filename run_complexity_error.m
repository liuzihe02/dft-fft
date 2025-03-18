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
num_expt = 10;  % Change this number as needed

% Preallocate matrices to store runtime measurements and RMSE values.
% Rows correspond to different sequence lengths; columns correspond to experiment repetitions.
runtime_my_dft_all = zeros(length(N_values), num_expt);
runtime_my_fft_all = zeros(length(N_values), num_expt);
runtime_mat_fft_all = zeros(length(N_values), num_expt); % Added for MATLAB's built-in FFT
rmse_my_dft_all    = zeros(length(N_values), num_expt);
rmse_my_fft_all    = zeros(length(N_values), num_expt);
rmse_mat_fft_all = zeros(length(N_values), num_expt);    % Added for MATLAB's built-in FFT

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
        runtime_my_dft_all(idx, expt) = toc;
        
        % Measure runtime for the custom fft_vectorized function
        tic;
        fft_signal = fft_vectorized(signal);
        runtime_my_fft_all(idx, expt) = toc;

        % Measure runtime for MATLAB's built-in FFT function
        tic;
        mat_fft_signal = fft(signal);
        runtime_mat_fft_all(idx, expt) = toc;

        % --- Compute reconstruction using MATLAB's IDFT (ifft) ---
        rec_my_dft = ifft(dft_signal);
        rec_my_fft = ifft(fft_signal);
        rec_mat_fft = ifft(mat_fft_signal);

        % --- Compute RMSE between original signal and reconstructed signal ---
        % Since ifft gives us complex values, we need to take the absolute value of the difference first before doing rmse
        % the . dot operator tells matlab to perform element wise operations
        rmse_my_dft_all(idx, expt) = sqrt(mean(abs(signal - rec_my_dft).^2));
        rmse_my_fft_all(idx, expt) = sqrt(mean(abs(signal - rec_my_fft).^2));
        rmse_mat_fft_all(idx, expt) = sqrt(mean(abs(signal - rec_mat_fft).^2));
    end
end

% Compute mean and standard deviation for each sequence length (across runs)
mean_runtime_my_dft = mean(runtime_my_dft_all, 2); %calculate mean across rows
std_runtime_my_dft  = std(runtime_my_dft_all, 0, 2); %std across rows
mean_runtime_my_fft = mean(runtime_my_fft_all, 2);
std_runtime_my_fft  = std(runtime_my_fft_all, 0, 2);
mean_runtime_mat_fft = mean(runtime_mat_fft_all, 2);
std_runtime_mat_fft  = std(runtime_mat_fft_all, 0, 2);

% Compute mean and standard deviation of RMSE across experiments for each sequence length
mean_rmse_my_dft = mean(rmse_my_dft_all, 2);
std_rmse_my_dft  = std(rmse_my_dft_all, 0, 2);
mean_rmse_my_fft = mean(rmse_my_fft_all, 2);
std_rmse_my_fft  = std(rmse_my_fft_all, 0, 2);
mean_rmse_mat_fft = mean(rmse_mat_fft_all, 2);
std_rmse_mat_fft  = std(rmse_mat_fft_all, 0, 2);

% %% PLOTTING THE COMPLEXITY RESULTS WITH ERROR BARS
% fig1 = figure;
% errorbar(N_values, mean_runtime_my_dft, std_runtime_my_dft, '-o', 'LineWidth', 2);
% hold on;
% errorbar(N_values, mean_runtime_my_fft, std_runtime_my_fft, '-s', 'LineWidth', 2);
% errorbar(N_values, mean_runtime_mat_fft, std_runtime_mat_fft, '-^', 'LineWidth', 2);
% xlabel('Sequence Length (N)');
% ylabel('Runtime (seconds)');
% title('Algorithmic Complexity of DFT and FFT');
% legend('my\_DFT', 'my\_FFT', 'MATLAB\_FFT', 'Location', 'NorthWest');
% grid on;
% set(gca, 'XScale', 'log', 'YScale', 'log');
% hold off;
% drawnow;  % ensure the figure is fully rendered
% print(fig1, 'complexity.png', '-dpng', '-r300', '-painters');
% fprintf('Algorithm complexity plot saved\n');

%% COMPLEXITY ANALYSIS
% Prepare data for complexity analysis
mean_runtime_data = {mean_runtime_my_dft, mean_runtime_my_fft, mean_runtime_mat_fft};
std_runtime_data = {std_runtime_my_dft, std_runtime_my_fft, std_runtime_mat_fft};
alg_names = {'my-DFT', 'my-FFT', 'MATLAB-FFT'};

% Analyze complexity
results = fit_complexity(N_values, mean_runtime_data, alg_names);

% Plot complexity results
plot_complexity(N_values, mean_runtime_data, std_runtime_data, alg_names, results);


%% PLOTTING RMSE (RECONSTRUCTION ERROR) WITH ERROR BARS
fig2 = figure;
plot_colors = get(gca, 'ColorOrder'); % Get the same default color order
markers = {'o', 's', '^'};
errorbar(N_values, mean_rmse_my_dft, std_rmse_my_dft, '-o', 'LineWidth', 2, 'Color', plot_colors(1,:));
hold on;
errorbar(N_values, mean_rmse_my_fft, std_rmse_my_fft, '-s', 'LineWidth', 2,'Color', plot_colors(2,:));
errorbar(N_values, mean_rmse_mat_fft, std_rmse_mat_fft, '-^', 'LineWidth', 2, 'Color', plot_colors(3,:));
xlabel('Sequence Length (N)');
ylabel('RMSE');
title('Reconstruction Error (RMSE)');
legend('my-DFT reconstruction', 'my-FFT reconstruction', 'MATLAB-FFT reconstruction', 'Location', 'NorthWest');
grid on;
set(gca, 'XScale', 'log');
hold off;
drawnow;  % ensure rendering
print(fig2, 'rmse.png', '-dpng', '-r300', '-painters');
fprintf('RMSE analysis plot saved\n');
