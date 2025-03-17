% Load the combined sine wave data for the corresponding input
% the variables (data_ch1,data_ch2) will be automatically loaded
% you can use the command 'whos' to check what variables are currently in the workspace
load('data/sine5_F0A1_300Hz.mat');

% Create time vector for plotting
%real_rate should be set at the sampling frequency
real_rate = 300;
%create a time vector starting from 0, step of 1/real_rate, all the way up till ending time
%ending time is determined by looking at .mat data
% make sure this t vector matches up with data_ch1
t = [0:1/real_rate:(length(data_ch1)-1)/real_rate]';

%set the input; make sure this aligns with the output data results
input_freq=5;
%this is the input we send to our force transducer
data_out = sin(2*pi*input_freq*t);

%we only need to analyse channel 2 since channel 1 is the input force transducer
f_ch2 = fft_vectorized(data_ch2);

%number of data points, because data may be padded with extra zeros
N_fft = length(f_ch2);
% array of 0 till number of data points
f = (0:N_fft-1)' * (real_rate / N_fft);


figure
hold on
plot(f, abs(f_ch2));

title('Transform of channel 2')
% V = axis;
% axis([0 50 V(3) V(4)])

% Save the most recent figure with explicit renderer settings, if some renders dont work, try either opengl or painters
%print('-dpng', '-r300', '-opengl', 'myanalyse.png')
print('-dpng','-r300','-painters','myanalyse.png');

fprintf('Plot of analysis saved in the current directory\n')

hold off



function X = dft_vectorized(x)
    % Compute the Discrete Fourier Transform (DFT) of input vector x using vectorized operations
    % x : input signal (vector)
    % X : DFT of x

    % Convert input to a column vector for consistency
    x = x(:);
    N = length(x);       % Number of samples in the signal
    
    % Create index vectors: n as a row vector and k as a column vector
    n = 0:N-1;           % Time indices (row vector)
    k = n';              % Frequency indices (column vector)
    
    % Construct the DFT matrix using the formula: exp(-1j*2*pi*k*n/N)
    W = exp(-1j * 2 * pi * k * n / N);
    
    % Multiply the DFT matrix with the signal to obtain the transform
    X = W * x;
end

function X = dft_loop(x)
    % Compute the Discrete Fourier Transform (DFT) of input vector x using naive loops
    % x : input signal (vector)
    % X : DFT of x

    % Convert input to a column vector for consistency
    x = x(:);
    N = length(x);       % Number of samples in the signal
    
    % Pre-allocate the output vector for efficiency
    X = zeros(N, 1);
    
    % Loop over each frequency bin k (from 0 to N-1)
    for k = 0:N-1
        sum_val = 0;
        % Loop over each time sample n (from 0 to N-1)
        for n = 0:N-1
            % Compute and accumulate the contribution for the k-th frequency component
            sum_val = sum_val + x(n+1) * exp(-1j * 2 * pi * k * n / N);
        end
        % MATLAB uses one-based indexing, so assign to X(k+1)
        X(k+1) = sum_val;
    end
end

function X = fft_vectorized(x)
    %  computes the radix-2 FFT recursively using vectorized operations.
    %   If the length of x is not a power of 2, it is zero-padded to the next power of 2, added at the end
    
    x = x(:);         % Ensure x is a column vector
    N = length(x);
    
    % If N is not a power of 2, zero-pad x to the next power of 2
    M = 2^nextpow2(N);  % nextpow2 returns the exponent so that 2^exponent >= N
    if M ~= N
        %add zeros to the end
        x = [x; zeros(M - N, 1)];
        N = M;  % Update N to the new length
    end
    
    % Base case: if the input length is 1, return x
    if N == 1
        X = x;
        return;
    end
    
    % Recursively compute FFT for even and odd indices
    X_even = fft_vectorized(x(1:2:end)); %select all the even indices
    X_odd  = fft_vectorized(x(2:2:end)); % selects all the odd indices
    
    % Compute twiddle factors (complex exponentials) in a vectorized manner
    factor = exp(-1j * 2 * pi * (0:N/2-1).' / N);
    
    % Combine the FFTs of the even and odd parts using the butterfly operation
    X = [X_even + factor .* X_odd;
         X_even - factor .* X_odd];
end
