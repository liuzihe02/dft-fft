function [ signal ] = sine( varargin )
%SINE Function to produce a single sine wave, or a combination of sine waves with different frequencies.
%   Produces a sine wave of frequency 'freq' with a resolution defined
%   by 't' (set in 'setiorates')
%
% E.g. data_out=sine(5,10); % this produces a signal with two componenets, at 5Hz and 10 Hz.  They are in phase.

global t;
signal = 0;

for k = 1:length(varargin)
    freq(k) = varargin{k}(1); % Cell array indexing
    signal = signal + sin(2*pi*t*freq(k));
end

signal = signal/length(varargin);
