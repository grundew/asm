%% Plane wave normal incidence frequency response

%% Chirp pulse
fs = 50e6;
tend = 50e-6;
t = (0:1/fs:tend)';
f0 = 50e3;
f1 = 600e3;

alpha = (f1-f0)/tend;
wndw = hann(length(t));
%y = wndw.*1j.*exp(1j*2*pi*alpha*t.^2/2);
y = chirp(t, f0, t(end), f1, 'linear', 270);

% Filter to remove the DC-component
n = 100;
wn = [10e3, 1000e3]./fs/2;
b = fir1(n, wn);
y = filtfilt(b, 1, y);

%% TODO
nfft = 2^18;
f = fftshift((-nfft/2:nfft/2-1)*fs/nfft);
Y = fft(y, nfft);

plotIR(t, y);

%% Fluid-fluid-fluid model
fluid1 = struct('v', 1500, 'density', 1000);
fluid3 = fluid1;
layer = struct('v', 5900, 'density', 7850);
theta = 0;
d = 12.4e-3;
[V, W] = fluidLayerReflectionCoefficient(f, theta, fluid3, layer, fluid1, d);

%% Plot reflection and transmission coefficient
ax = plotIR(f, V);
plotIR(f, W, true, ax);
legend(ax(1), 'Reflection', 'Transmission')

%% Compare pulse spectrum and reflection coefficient
figure;
plotyy(f, abs(Y), f, abs(V));
legend('Reflection', 'FFT(pulse)')

%% Convolve the pulse and the reflection coefficient
infft = length(t);
x = ifft(Y.*V, infft);
plotIR(t, x, true);

%% Impulse response of the reflection coefficient
v = ifft(V, nfft);
ax = plotIR(0:length(v)-1, v);