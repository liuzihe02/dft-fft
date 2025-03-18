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
    % row array
    factor = exp(-1j * 2 * pi * (0:N/2-1).' / N);
    
    % Combine the FFTs of the even and odd parts using the butterfly operation
    X = [X_even + factor .* X_odd;
         X_even - factor .* X_odd];
end