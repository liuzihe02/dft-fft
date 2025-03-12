# dft-fft
Investigating the Discrete Fourier Transform and Fast Fourier Transform, with applications to response of buildings in earthquakes

## Running A4 Lab Remotely

I have set up the A4 lab `matlab` code in the `/a4` folder and relevant data in `/data`. To generate figures as in the A4 lab, run the following matlab code (assuming we are in the home directory):

```matlab
% Load the combined sine wave data
load('sine_5_10Hz.mat');

% Create time vector
t = [0:1/real_rate:(length(data_ch1)-1)/real_rate]';

% Create output signal for reference
data_out = sin(2*pi*5*t) + sin(2*pi*10*t);

% Display signals
analyse1;

% Save figure
saveas(gcf, 'fig1c.pdf');
```
