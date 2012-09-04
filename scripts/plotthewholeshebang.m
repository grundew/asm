%% Plot the chirp with the spectrum
figure
subplot(211)
plot(t, real(y))
subplot(212)
plot(f, db(abs(Y)/max(abs(Y))), '.')

%% Plot the frequency and impulse response of the observation point
figure
subplot(211)
plot(f, real(p), '.')
title('Frequency response of the plate')
subplot(212)
plot(f, imag(p), '.')

figure
plot((0:nfft-1)/fs, real(fft(p, nfft)))
title('Impulse response of the plate')

%% Plot the convolved pulse
figure
plot(tt, real(yy))
xlabel('time')
ylabel('amplitude')
title('Rx pulse (through transmission)')