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
f_ch2 = fft(data_ch2);
x=t*real_rate/(t(length(t)));

figure
hold on

plot(x,abs(f_ch2));
title('Transform of channel 2')
% V = axis;
% axis([0 50 V(3) V(4)])

% Save the most recent figure with explicit renderer settings, if some renders dont work, try either opengl or painters
%print('-dpng', '-r300', '-opengl', 'myanalyse.png')
print('-dpng','-r300','-painters','myanalyse.png');

fprintf('Plot of analysis saved in the current directory\n')

hold off



% %% Custom DFT Function
% % This function computes the Discrete Fourier Transform (DFT) of an input vector.
% function X = my_dft_vectorized(x)
%     N = length(x);       % Number of samples in the signal
%     n = 0:N-1;           % create a row array of Time indices
%     % ' is the transpose operation
%     k = n';              % Frequency indices as a column vector
%     % Create the DFT matrix using the formula exp(-j*2*pi*k*n/N)
%     % this is an outer product, which gives the DFT matrix
%     W = exp(-1j * 2 * pi * k * n / N);
%     X = W * x(:); % Multiply the DFT matrix by the input signal (converted to a column vector), the (:) operator converts it to a column vector
% end
