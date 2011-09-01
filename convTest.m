clear classes
%% Generate a complex chirp and its fft
fs = 50e6;
tend = 50e-6;
t = (0:1/fs:tend)';
f0 = 50e3;
f1 = 600e3;

alpha = (f1-f0)/tend;
wndw = hann(length(t));
y = wndw.*1j.*exp(1j*2*pi*alpha*t.^2/2);

nfft = 2^16;
f = (-nfft/2:nfft/2-1)*fs/nfft;
f = fftshift(f);
Y = fft(y, nfft);

% Pick out the necessary frequencies. All others are outside the bandwidth
% of the pulse (hence = zero).
fmin = 0;
fmax = 1e6;
id = find(f<fmax & f>fmin);
ff = f(id);

%% Do the whole she bang
q = linspace(0, sin(50*pi/180), 200);
xo = [0 0.1];
z = 10e-2;
d = 12.4e-3;
a = 3e-3;

wni = WaveNumberIntegration(xo, z, q, a, MultiLayerModel('watersteelwater', d));
[Pr, Pt] = wni.doAll(ff);

%% Plotting & comparing with normal incidence
figure
ax1 = subplot(211);
ax2 = subplot(212);
hold(ax1, 'all')
hold(ax2, 'all')

yPulse = ifft(Y(id).*Pr(1, :)', nfft);
tailPrd = (abs(fft(yPulse(2400:end), nfft)).^2)/nfft;

plot(ax1, real(yPulse)/max(real(yPulse(:))))
plot(ax2, f, tailPrd/max(tailPrd(:)))
xlim([0 500e3]);

[R, T] = waterSteelWaterC1(ff, 0, d);
R(f==0) = 0;
yPulse_norm = ifft(Y(id).*R, nfft);

plot(ax1, real(yPulse_norm)./max(real(yPulse_norm(:))))
tailPrd = (abs(fft(yPulse_norm(2400:end), nfft)).^2)/nfft;
plot(ax2, f, tailPrd/max(tailPrd(:)))
legend(ax1, 'Finite source', 'Plane wave, normal incidence')