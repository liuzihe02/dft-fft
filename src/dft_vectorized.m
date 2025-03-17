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
