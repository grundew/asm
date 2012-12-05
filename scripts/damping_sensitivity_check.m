%% Test the use of damping factor in the speed of sound

%% Generate the chirp pulse
fs = 2e6;
tend = 50e-6;
t = (0:1/fs:tend)';
f0 = 200e3;
f1 = 800e3;

% Real chirp
wndw = gausswin(length(t));
y = wndw.*chirp(t, f0, tend, f1, 'linear', 270);
y = [zeros(100, 1); y; zeros(100, 1)];
t = (0:length(y)-1)/fs';
y = conj(hilbert(y));
y = y - mean(y);

nfft = 2^18;
f = (0:nfft-1)*fs/nfft;
Y = ifft(y, nfft);

%% Do the whole she bang
rho_fluid = 1000;
v_fluid = 1500;
v = 0;
v_layer = 5850;
theta = 0.1;

%% Fluid-fluid-fluid model
fluid1 = struct('v', v_fluid, 'density', rho_fluid);
fluid3 = fluid1;
layer = struct('v', v_layer, 'density', 7850, 'vShear', 3218);

d = 10e-3;

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, d);
[~, W] = analyticRTFast(f, theta, model);
x = fft(Y(:).*W(:), nfft);
xprd = abs(ifft(x, nfft)).^2;
tailstart = 800;
tailend = 1200;
xtailprd = abs(ifft(x(tailstart:tailend), nfft)).^2;
tx = (0:nfft-1)/fs;

%% Plot the impulse response of the transmission and reflection coefficients
figure
subplot(211)
plot((0:nfft-1)/fs, abs(fft(W, nfft)));
xlabel('Time')
ylabel('Abs')
title('Impulse response')
subplot(212)
plot((0:nfft-1)/fs, unwrap(angle((fft(W, nfft)))))
ylabel('Angle')
xlabel('Time')

%% Plot tail signal and the whole signal
plotSignal(tx, x, tx(tailstart:end), x(tailstart:end),...
    f, xprd, xtailprd, 0.5/d*[(1:2)*layer.v, 3*layer.vShear]);