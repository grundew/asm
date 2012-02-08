%% Code for convolving a plane wave with normal incidence with a chirp
% 1) Fluid - fluid - fluid
% 2) Fluid - solid - fluid (at very small angle)

clear;
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

% Plot the pulse
plotIR(t, y, false);
% Plot the fft of the pulse
plotIR(f, Y, true);

%% Reflection and transimission coefficients: Fluid - solid - fluid model
qs = 0.001;
fluid.v = 1500;
fluid.density = 1000;
solid.v = 5900;
solid.vShear = 3150;
solid.density = 7850;
d = 12.4e-3;
model = MultiLayerModel(fluid, solid, fluid, d);
[Rs, Ts] = fluidSolidFluidReflectionCoefficient(ff, asin(qs), model);
ax = plotIR(ff, Rs, true);
plotIR(ff, Ts, true, ax);
legend(ax(1), 'R_fluid', 'T_fluid')

%% Reflection and transmission coefficients: Fluid - fluid - fluid model
qf = 0;
[Rf, Tf] = fluidLayerReflectionCoefficient(f, asin(qf), fluid, solid, fluid, d);

plotIR(f, Rf, true, ax);
plotIR(f, Tf, true, ax);
legend(ax(1), 'R_fluid', 'T_fluid', 'R_solid', 'T_solid')
%% Plotting & comparing with normal incidence
[yPulses, ts] = convTimeFreq(y, fs, nfft, Rs, ff, true);
[yPulsef, tf] = convTimeFreq(y, fs, nfft, Rf, f, true);

axt = plotIR(tf, yPulsef, false);
plotIR(ts, yPulses, false, axt);

% yPulse = ifft(Y(id).*Tf, nfft);
% tailPrd = (abs(fft(yPulse(2400:end), nfft)).^2)/nfft;
% 
% plot(ax1, real(yPulse)/max(real(yPulse(:))))
% plot(ax2, f, tailPrd/max(tailPrd(:)))
% xlim([0 500e3]);
% 
% yPulse_norm = ifft(Y(id).*Rf, nfft);
% 
% plot(ax1, real(yPulse_norm)./max(real(yPulse_norm(:))))
% tailPrd = (abs(fft(yPulse_norm(2400:end), nfft)).^2)/nfft;
% plot(ax2, f, tailPrd/max(tailPrd(:)))
% legend(ax1, 'Solid layer', 'Fluid layer')