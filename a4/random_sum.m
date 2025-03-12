function [ signal ] = random_sum( )
%RANDOM_SUM Function to generate a superposition of sine waves with 'random'
%phases
%   Produces a superposition of sine waves with a range of frequencies, with 'random' phases, and with a resolution defined
%   by 't' (set in 'setiorates')

global t;
global real_time;
global real_rate;

x=0;

for i =1/real_time:1/real_rate:real_rate/2;
    phi=2*pi*rand;
    x = x+sin(2*pi*i*t+phi);
end

x = 10*x./ ( (real_rate/2 - 1/real_time) / (1/real_rate));

signal = x;
