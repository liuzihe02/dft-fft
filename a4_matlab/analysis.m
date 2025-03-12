% analysis1 : a simple analysis of the data, takes ffts and plots the
% transfer function (can do this for either A/F, or A over the output
% signal from the card - results seem very similar, as expected)

% works well with a data output of chirp(t,0,10,20), for t = 0:0.01:10, but
% transpose it!!

%Plot the time domain response of data sent to the card and recieved from
%the card
plot(t,force_signal);
hold on
plot(t,data_out);

%Plot the frequency response of data sent to the card and recieved from the
%card
figure
Fourier_Transform(force_signal,acceleration_signal,data_out);
