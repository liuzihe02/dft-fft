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