%% Plane wave normal incidence frequency response

%% Chirp pulse
fs = 2e6;
tend = 50e-6;
t = (0:1/fs:tend)';
f0 = 50e3;
f1 = 800e3;

alpha = (f1-f0)/tend;
wndw = hann(length(t));
y = wndw.*exp(-1j*2*pi*alpha*t.^2/2);

n = 5;
[z, p, k] = butter(n, 2*f0/fs, 'high');
[sos, g] = zp2sos(z,p,k);
y = filtfilt(sos, g, y);

nfft = 2^10;
f = fftshift((-nfft/2:nfft/2-1)*fs/nfft);
idfneg = f<0;
idfpos = f>0;
Y = ifft(y, nfft);

%% Fluid-fluid-fluid model
% fluid1 = struct('v', 3500, 'density', 1000);
fluid1 = struct('v', 340, 'density', 1.2);
fluid3 = fluid1;
layer = struct('v', 5850, 'density', 7850, 'vShear', 3218);
theta = 1e-2;
d = 10.15e-3;

% [V, W] = fluidLayerReflectionCoefficient(f, theta, fluid3, layer, fluid1, d);

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, d);
V = fluidSolidFluidReflectionCoefficient(f, theta, model);

%% Convolve the pulse and the reflection coefficient
% V = zeros(size(Y));
% Y(idfneg) = 0;
% V(idfpos) = R;
% V(idfneg) = fliplr(conj(R));
% V(f==0) = 0;

V = V(:);
Y = Y(:);
x = fft(Y.*V, nfft);
xprd = abs(ifft(x, nfft)).^2;
xtailprd = abs(ifft(x(300:end), nfft)).^2;
tx = (0:nfft-1)/fs;
figure
subplot(211)
plot(tx, real(x));
subplot(212)
plot(f(idfpos), db(xprd(idfpos), 'P'))
xlim([0 800e3])
figure
plot(f(idfpos), xtailprd(idfpos))
% %% Plot reflection and transmission coefficient
% ax = plotIR(f, V);
% plotIR(f, W, true, ax);
% legend(ax(1), 'Reflection', 'Transmission')
% 
% %% Compare pulse spectrum and reflection coefficient
% figure;
% plotyy(f, abs(Y), f, abs(V));
% legend('Reflection', 'FFT(pulse)')
% 
% %% Impulse response of the reflection coefficient
% v = fft(V, nfft);
% ax = plotIR(0:length(v)-1, v);