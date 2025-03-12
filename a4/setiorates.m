function [] = setiorates(time,rate)
% setiorates - sets the input and output sampling rates and time.
% Call the function with the following arguments:
% time, in seconds, to take samples for, on the input channels;
% desired sampling rate, in samples per second (this is used for input and
% output)
% NB: the time for the output signal is determined by the length of the
% data_out variable (data points/sample rate).
%
% E.g. setiorates(10,1000) % 10 seconds, 1kHz sampling rate

% only works IN PERSON

global in;
global out;
global real_time;
global real_rate;
global t;

if ~exist('time') || ~exist('rate')
    msgbox ('Please supply a desired sampling time and sampling rate','Function argument problem','warn')
    return
end
if (time < 0.1) || (time > 600) || (rate < 0.6) || (rate > 250000)
    msgbox ('Please enter a valid time and sampling rate.  Valid times are from 0.1 to 600 seconds (the first function argument), valid sample rates are from 0.6 Hz to 250 kHz (the second function argument).','Function argument problem','warn')
    return
end

%% Input setup

set(in,'SamplesPerTrigger',time*rate)
set(in,'SampleRate',rate)

%% Output setup

set(out,'SampleRate',rate)

%% Set up a vector 't' for use in the chirp function (useful for generating data_out)

datalength=get(in,'SamplesPerTrigger');
datarate=get(in,'SampleRate');

real_time=datalength/datarate;
real_rate=datarate;

t=0:1/datarate:(datalength/datarate-1/datarate);
t=t';
