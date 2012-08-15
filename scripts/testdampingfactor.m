%% Test the use of damping factor in the speed of sound

%% Generate the chirp pulse
fs = 2e6;
tend = 50e-6;
t = (0:1/fs:tend)';
f0 = 150e3;
f1 = 350e3;

% Real chirp
wndw = gausswin(length(t));
y = wndw.*chirp(t, f0, tend, f1, 'linear', 270);
y = [zeros(100, 1); y; zeros(100, 1)];
t = (0:length(y)-1)/fs';
% Analytic chirp
% alpha = (f1-f0)/tend;
% wndw = gausswin(length(t));
% y = wndw.*exp(-1j*2*pi*alpha*t.^2/2);

% Filter the pulse
% n = 10;
% [z, p, k] = butter(n, 2*f0/fs, 'high');
% [sos, g] = zp2sos(z,p,k);
% y = filtfilt(sos, g, y);

nfft = 2^18;
f = fftshift((-nfft/2:nfft/2-1)*fs/nfft);
Y = ifft(y, nfft);

%% Do the whole she bang
rho_fluid = 1000;
v_fluid = 1500;
v = 0;
v_layer = 5850 - 500i;
theta = 0;

%% Fluid-fluid-fluid model
fluid1 = struct('v', v_fluid, 'density', rho_fluid);
fluid3 = fluid1;
layer = struct('v', v_layer, 'density', 7850, 'vShear', 3218);
%theta = 1e-4;
d = 10.15e-3;

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, d);

tic
[V, W] = analyticRT(f, theta, model);
toc
x = fft(Y(:).*V(:), nfft);
xprd = abs(ifft(x.*hann(length(x)), nfft)).^2;
tailstart = 800;
tailend = 1200;
xtailprd = abs(ifft(x(tailstart:tailend).*hann(tailend-tailstart+1), nfft)).^2;
tx = (0:nfft-1)/fs;

%% Plot the impulse response of the transmission and reflection coefficients
figure
subplot(211)
plot((0:nfft-1)/fs, real(fft(V, nfft)));
subplot(212)
plot((0:nfft-1)/fs, real(fft(W, nfft)));

%% Plot tail signal and the whole signal
plotSignal(tx, x, tx(tailstart:end), x(tailstart:end),...
    f, xprd, xtailprd, 0.5/d*[(1:2)*layer.v, 3*layer.vShear]);