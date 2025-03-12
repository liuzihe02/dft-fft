% This script runs the data capture board in input-output mode, using the data in data_out.

%% checks data_out and adds a zero on the end (to make the board output zero after it's been used)
data_out(length(data_out))=0;

if length(data_out)==get(in,'SamplesPerTrigger') % checks that data_out has the right length

% A run to set the equipment up before recording the data that will be
% used.  This ensures that, when the data is logged, the signal is close
% to periodic.

%% runs the board and gets the data
samples=get(in,'SamplesPerTrigger');
rate=get(in,'SampleRate');
time=samples/rate;
setiorates(time*2,rate);

putdata(out,[data_out;data_out])
start([in out])
trigger([in out])
data_in = getdata(in);

data_ch1=data_in(:,1);
data_ch2=data_in(:,2);

% discard the first half of the data
l=length(data_out);
data_ch1=data_ch1(l+1:l*2);
data_ch2=data_ch2(l+1:l*2);

data_ch1=data_ch1-mean(data_ch1);
data_ch2=data_ch2-mean(data_ch2);

setiorates(time,rate);

else
    display('Please create data_out with a vaild set of run_time and sampling rate.  You probably need to set data_out again after changing the time and/or sampling rate.');
end
