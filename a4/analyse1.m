figure
hold on

plot(t,data_out,t,data_ch1,t,data_ch2)
title('Time-domain plot')
xlabel('Time (seconds)')
ylabel('Signal (volts)')
legend('Generated Signal','Channel 1','Channel 2')

hold off
