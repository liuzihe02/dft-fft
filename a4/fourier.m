function [ transform ] = fourier( signal )
%FOURIER Summary of this function goes here
%   Detailed explanation goes here

%s1 = cputime;

transform = fft(signal);

% m=1:1:length(signal);
% k=1:1:length(signal);


%Time_taken = cputime - s1
