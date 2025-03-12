function [ signal ] = sweep( f1, f2 )
%SWEEP Function to generate a sine sweep
%   Produces a sine sweep between frequencies f1 and f2 with a resolution defined
%   by 't' (set in 'setiorates')
global t;

signal = chirp(t,f1,t(length(t)),f2);