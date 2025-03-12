% initialise - initialises global variables and sets up the board for simultaneous input and output

clear all

global in out real_time real_rate t;
real_time = 0;
real_rate=0;
t=0;

%% Input setup
names = {'channel 1';'channel 2'};

in = analoginput('nidaq','Dev1');
addchannel(in,0:1,names)
set(in,'InputType','SingleEnded')
set(in,'TriggerType','manual')

%% Output setup
names = {'output 1 (channel 9)'};

out = analogoutput('nidaq','Dev1');
addchannel(out,0,names)
set(out,'TriggerType','manual')

%% prepare board and data

set([in out],'StopFcn',@daqcallback)
