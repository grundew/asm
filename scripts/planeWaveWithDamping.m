% This script demonstrates the use of planeWaveTransmissionTimeSignal and
% planeWaveTimeSignal to simulate a single plane wave with a given angle of
% incidence. The transmission coefficient for the solid plate is calculated
% and convolved with the excitition pulse to give the time signal.

%% Pulse
fs = 2e6;
nfft = 2^21;
tend = 50e-6;
t = (0:1/fs:tend)';
f0 = 200e3;
f1 = 800e3;
wndw = gausswin(length(t));
y = wndw.*chirp(t, f0, tend, f1, 'linear', 270);
y = conj(hilbert(y));

%% Material properties
alphaLambda = 10^(0.008/20);
% alphaLambda = 1;
fluid1 = struct('v', 350, 'density', 1.5);
layer = struct('v', 5850, 'density', 7850, 'vShear', 3218);
model = MultiLayerModel(fluid1, layer, fluid1, 10e-3);
dist = 10e-2;
theta = 0;

%% Compute transmitted pulse #1
tic
[xT, t, T] = planeWaveTransmissionTimeSignal(model, y, t, theta, dist, alphaLambda, nfft);
toc

%% Compute transmitted pulse #2
tic
[~, xT2, t2, ~, T2] = planeWaveTimeSignal(model, y, t, theta, dist, nfft);
toc

%% Plot time signals
figure
plot(t, real(xT), '.', t2, real(xT2), 'o')
title('Transmission')
legend('With damping', 'Without damping')
xlabel('Time')
ylabel('Amplitude')

%% Plot frequency spectrum
f = (0:nfft-1)*fs/nfft;
fres = 0.5*5850/10e-3;
figure
plot(f/fres, db(1/nfft*abs(ifft(xT, nfft))), '.', f/fres, db(1/nfft*abs(ifft(xT2, nfft))), 'o')
xlabel('Normalized frequency')
ylabel('PSD')
title('Frequency spectrum')

%% Plot the transmission coefficient
figure
plot(f/fres, abs(T), '.', f/fres, abs(T2), 'o')
xlabel('Frequency (Normalized to resonance)')
ylabel('abs T')
title('Transmission coefficient')
legend('With damping', 'Without damping')

%% Plot the decay of the time signal
figure
semilogy(t, abs(xT), t, abs(xT2))
xlabel('Time')
ylabel('Abs amplitude')
title('Signal decay')
legend('With damping', 'Without damping')
