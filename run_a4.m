%% IMPORT STATEMENTS
% Add the functions folder `src` to the MATLAB search path
addpath(fullfile('src'));

%% LOAD AND PROCESS DATA

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

%% DFT/FFT TRANSOROMATION

%we only need to analyse channel 2 since channel 1 is the input force transducer
f_ch2 = fft_vectorized(data_ch2);

%% PLOT DATA

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
