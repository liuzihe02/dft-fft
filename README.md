# Investigation of DFT and FFT
Investigating the Discrete Fourier Transform (DFT) and Fast Fourier Transform (FFT), with applications to response of buildings in earthquakes.

## A4 Lab

We run the A4 lab remotely to speed up gathering and processing of data. `F0` is the force transducer on the oscillator, `A1;A2;A3` are the sensors on floors 1,2,3 of the structure, respectively. `data_ch1` is the data from the input (force transducer), while `data_ch2` is the reading from the accelerometer of the respective floor.

### Data Format

1. files of the form  `sineX_F0AZ_N.mat`  : which give the data resulting from  outputting a sine wave of frequency X Hz to the card (and therefore to  the structure) and  then interrogating the response at points F0 (the force transducer) and the sensor AZ, where Z=1,2,3 and represents the number of the floor. The .mat file contains time series data from these two channels (channel 1 is F0, channel 2 is AZ). Each file contains 15seconds of data. See Lab Sheet B for more detail.
2. files of the form  `sineX_Y_F0AZ_N.mat`  : which give the data resulting from  outputting a sine wave of frequency X Hz plus a sine wave of frequency Y to  the card (otherwise as in 1.).
3. files of the form `sweep1to20_F0AZ_N.mat`  : which give the data resulting from  outputting a chirp signal between 1 and 20 Hz (otherwise as in 1.). See Lab Sheet B for more detail.
4. files of the form `random20timesLarger_F0AZ_N.mat`  : which give the data resulting from  outputting a superposition of sine waves with random phases (otherwise as in 1. but data has been amplified by a factor of 20).

## Repo Structure

Our implementations of DFT and FFT are stored in the `src` folder. Scripts titled `run_a4.m` call these functions for data processing and plotting accordingly.

### Running Experiment A4

I have set up the A4 lab `matlab` code in the `/a4` folder and relevant data in `/data`. The following contains a summarized script to run the A4 experiments remotely: (assuming we are in the home directory):

```matlab

% Load the combined sine wave data for the corresponding input
% the variables (data_ch1,data_ch2) will be automatically loaded
% you can use the command 'whos' to check what variables are currently in the workspace
load('data/sine5_F0A1_300Hz.mat');

% Create time vector for plotting
%real_rate should be set at the samplimg frequency
real_rate = 300;
%create a time vector starting from 0, step of 1/real_rate, all the way up till ending time
%ending time is determined by looking at .mat data
% make sure this t vector matches up with data_ch1
t = [0:1/real_rate:(length(data_ch1)-1)/real_rate]';

%set the input; make sure this aligns with the output data results
input_freq=5;
%this is the input we send to our force transducer
data_out = sin(2*pi*input_freq*t);

% run our analysis
run('a4/analyse1.m')
% Save the most recent figure with explicit renderer settings, using opengl
print('-dpng', '-r300', '-opengl', 'analyse1.png')
fprintf('Plot of analysis1 saved in the current directory\n')

% run our analysis
run('a4/analyse2.m')
% Save the most recent figure with explicit renderer settings, using opengl
print('-dpng', '-r300', '-opengl', 'analyse2.png')
fprintf('Plot of analysis2 saved in the current directory\n')
```

## Notes

- In `MATLAB`, all functions must come after executable scripts, at the end
- `MATLAB` uses 1-based indexing!
