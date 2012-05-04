%% Parameters
a = 12e-3;
h = 10e-2;
fluid.v = 350;
fluid.density = 1.4;
solid.v = 5850;
solid.vShear = 3218;
solid.density = 7850;
d = 10.15e-3;
model = MultiLayerModel(fluid, solid, fluid, d);

%% Test time signal
% Chirp pulse
fs = 2e6;
nfft = 2^12;
tend = 50e-6;
ty = (0:1/fs:tend)';
f0 = 50e3;
f1 = 800e3;

alpha = (f1-f0)/tend;
wndw = hann(length(ty));
y = wndw.*exp(-1j*2*pi*alpha*ty.^2/2);

% Filter to remove the DC-component
n = 5;
[z, p, k] = butter(n, 2*f0/fs, 'high');
[sos, g] = zp2sos(z,p,k);
y = filtfilt(sos, g, y);

%% Calculate time signal
nq = 2^15;
wni = WaveNumberIntegration2(a, h, [], model, nq);
% wni = WaveNumberIntegration2(a, h, [], model, 1, q);

nx = 1;
[x, tx] = calculateReflectedPressureRx(wni, a, nx, y, ty, nfft);

%% Plot the tail
tailstart = 400;
taillength = 600;
nfft = 2^13;
f = fftshift((-nfft/2:nfft/2-1)*fs/nfft);
Xprd = 1/nfft*abs(ifft(x, nfft)).^2;
xtail = x(tailstart:tailstart+taillength-1);
Xtailprd = 1/nfft*abs(ifft(xtail, nfft)).^2;

% Pick out positive frequencies
ixfpos = f>0;
f = f(ixfpos);
Xprd = Xprd(ixfpos);
Xtailprd = Xtailprd(ixfpos);

figure
ax(1) = subplot(211);
plot(f, Xprd)
ax(2) = subplot(212);
plot(f, Xtailprd)
linkaxes(ax, 'x')

%% Plot the time signal
figure
plot(tx, real(x))

%% Upsample
nup = 4*length(x);
y = interpft(x, nup);
dy = (tx(2)-tx(1))*length(x)/nup;
ty = 0:dy:(nup-1)*dy;

%% Save
save(sprintf('results_%s.mat', datestr(now, 'yyyymmddTHHMMSS')))
fprintf('Done, %s...\n', datestr(now, 'yyyymmddTHHMMSS'))