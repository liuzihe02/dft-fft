f_ch1 = fourier(data_ch1);
f_ch2 = fourier(data_ch2);
x=t*real_rate*1/t(length(t));

figure
hold on

plot(x,abs(f_ch2)./abs(f_ch1));
title('Transfer function')
xlim([0 x(length(x))/2])
xlabel('Frequency (Hz)')
ylabel('Response (ch2/ch1)')

hold off