%% Add functions folder 'src' to the MATLAB search path
addpath(fullfile('src'));

%% PARAMETERS
% Sampling frequency (Hz)
real_rate = 300;
% Input frequency for the ideal sine wave (Hz) - modify as needed
input_freq = 5;

% Define logarithmic range of sequence lengths (using powers of 2)
k_min = 5;    % smallest exponent (i.e., N = 2^5 = 32)
k_max = 14;   % largest exponent (i.e., N = 2^12 = 4096)
N_values = 2.^(k_min:k_max);  % sequence lengths

% Preallocate arrays to store runtime measurements (in seconds)
runtime_dft = zeros(size(N_values));
runtime_fft = zeros(size(N_values));

%% ALGORITHMIC COMPLEXITY ANALYSIS
% Loop over each sequence length
for idx = 1:length(N_values)
    N = N_values(idx);
    
    % Create a time vector for N samples
    t = (0:N-1)' / real_rate;
    
    % Generate the ideal sine wave signal
    signal = sin(2*pi*input_freq*t);
    
    % Measure runtime for the custom dft_vectorized function
    tic;
    dft_signal = dft_vectorized(signal);
    runtime_dft(idx) = toc;
    
    % Measure runtime for the custom fft_vectorized function
    tic;
    fft_signal = fft_vectorized(signal);
    runtime_fft(idx) = toc;
end

%% PLOTTING THE RESULTS
figure;
% Plot DFT runtime on a log-log scale
loglog(N_values, runtime_dft, '-o', 'LineWidth', 2);
hold on;
% Plot FFT runtime on a log-log scale
loglog(N_values, runtime_fft, '-s', 'LineWidth', 2);
xlabel('Sequence Length (N)');
ylabel('Runtime (seconds)');
title('Algorithmic Complexity: dft\_vectorized vs fft\_vectorized');
legend('dft\_vectorized', 'fft\_vectorized', 'Location', 'NorthWest');
grid on;
hold off;

% Save the plot as a PNG image with specified renderer settings
print('-dpng', '-r300', '-painters', 'complexity.png');
fprintf('Algorithm complexity plot saved as algorithm_complexity.png\n');
