%% Excitation pulse
nfft = 2^13;
f = (0:nfft-1)/nfft;
n = 512;
y = amgauss(n, n/2, n/2).*fmlin(n, 0, 0.5);
y = [zeros(64, 1); y; zeros(64, 1)];
n = length(y);
t = (0:n-1);
Y = ifft(y, nfft);

%% Plot the pulse
figure
subplot(211)
plot(t, real(y))
subplot(212)
plot(f, db(abs(Y).^2, 'power'))

%% Fluid-fluid-fluid model
theta = 0.01;
fluid1 = struct('v', 0.1, 'density', 0.1);
fluid3 = fluid1;
layer = struct('v', 1, 'density', 1, 'vShear', 0.55);
d = 2.5;
fres = 0.5/d;

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, d);

[R, T] = analyticRTFast(f, theta, model);

%% Convolve the pulse and the reflection coefficient
xt = fft(Y(:).*T(:), nfft);
xr = fft(Y(:).*R(:), nfft);
tx = 1:nfft;

%% Short time fourier transform
tfrstft(xt, tx, 2^12);