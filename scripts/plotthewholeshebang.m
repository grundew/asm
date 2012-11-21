%% Plot the chirp with the spectrum
figure
subplot(211)
plot(t, real(y))
subplot(212)
plot(f/fres, db(abs(Y)/max(abs(Y))), '.')

%% Extract the tail signal
tailstart = 1000;
taillength = 1000;
tailend = 700 + taillength;
ytail = real(yr(tailstart:tailend));
ttail = tt(tailstart:tailend);
Ytail = abs(ifft(ytail.*hann(length(ytail)), nfft)).^2/nfft;

%% Plot the tail signal and tail spectrum
figure
subplot(211)
plot(ttail, ytail)
xlabel('Time')
subplot(212)
plot(f/fres, db(Ytail)/max(Ytail), '.')
xlabel('f/f_{res}')

%% Plot the frequency and impulse response of the observation point
figure
subplot(211)
plot(f/fres, real(pt), '.')
title('Frequency response of the plate')
xlabel('f/f_{res}')
subplot(212)
plot(f/fres, imag(pt), '.')
xlabel('f/f_{res}')

figure
plot((0:nfft-1)/fs, real(fft(pt, nfft)))
title('Impulse response of the plate')
xlabel('Time')

%% Plot the convolved pulse
figure
plot(tt, real(yt))
xlabel('time')
ylabel('amplitude')
title('Rx pulse (through transmission)')