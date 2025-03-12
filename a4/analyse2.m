f_ch1 = fourier(data_ch1);
f_ch2 = fourier(data_ch2);
x=t*real_rate/(t(length(t)));

figure
hold on

subplot(1,2,1)
plot(x,abs(f_ch1));
title('Transform of channel 1')
V = axis;
axis([0 50 V(3) V(4)])

subplot(1,2,2)
plot(x,abs(f_ch2));
title('Transform of channel 2')
V = axis;
axis([0 50 V(3) V(4)])

hold off