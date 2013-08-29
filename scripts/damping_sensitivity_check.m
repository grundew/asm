% This is a script for investigating the transmission signal sensitivity to
% the damping factor.

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
% 0.008 dB/m
alphaLambda_dB = [0, 0.008, 0.08, 0.8, 1];
fluid1 = struct('v', 1500, 'density', 1000);
layer = struct('v', 5850, 'density', 7850, 'vShear', 3160);
model = MultiLayerModel(fluid1, layer, fluid1, 10e-3);
dist = 10e-2;
theta = 0;

f = (0:nfft-1)*fs/nfft;
fres = 0.5*5850/10e-3;
Y = ifft(y, nfft);
legendstr = {};

for i = 1:length(alphaLambda_dB)
    al = alphaLambda_dB(i);
    [xT, t, T] = planeWaveTransmissionTimeSignal(model, y, t, theta, dist, al, nfft);
    
    legendstr = cat(1, legendstr, num2str(al));
    
    %% Plot time signals
    figure(1)
    hold all
    plot(t, real(xT), '.')
    title('Transmission')
    xlabel('Time')
    ylabel('Amplitude')
    
    %% Plot frequency spectrum
    figure(2)
    hold all
    XT = 1/nfft*ifft(xT, nfft);
    plot(f/fres, db(abs(XT)))
    xlabel('Normalized frequency')
    ylabel('PSD')
    title('Frequency spectrum')
    
    %% Plot the transmission coefficient and pulse spectrum
    figure(3)
    hold all
    plot(f/fres, abs(T))
    xlabel('Frequency (Normalized to resonance)')
    ylabel('abs T')
    title('Transmission coefficient')
    
    %% Plot the decay of the time signal
    figure(4)
    hold all
    plot(t, db(abs(xT)))
    xlabel('Time')
    ylabel('Abs amplitude')
    title('Signal decay')
end

%% Add the excitation pulse spectrum and legends
figure(3)
hold all
plot(f/fres, abs(Y))
for i = 1:4
    figure(i)
    legend(legendstr)
end