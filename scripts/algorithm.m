% Algorithm

%% 1. Compute the chirp
fs = 15e6;
tend = 50e-6;
tChirp = 0:1/fs:tend;
f0 = 100e3;
fend = 500e3;
yChirp = hann(length(tChirp))'.*chirp(...
    tChirp, f0, tend, fend, 'linear', 270);

nfft = 2^11;
fx = fftshift((-nfft/2:nfft/2-1)*fs/nfft);

Y = fft(yChirp, nfft);
Ydb = 20*log10(abs(Y./max(abs(Y))));
% Plot everything above -65 dB
idp = find(Ydb > -65 & fx > 0);
idm = find(Ydb > -65 & fx < 0);

% Plot
plot(fx, Ydb)
hold all
plot([min(fx) max(fx)], [-65 -65])
plot(fx([idp, idm]), Ydb([idp, idm]), 'r.')

% Pick out positive frequencies
f = fx(idp);

%% 2. Calculate plane wave Fourier coefficients of the incoming pressure.
r = 0.12; % 12 cm
[X, kx, ky] = circularAperture(r, 0);
% Pick the centre line
Xc = X(ceil(size(X, 1)/2), :);

%% 3. Choose distance from transducers to plate, h1 & h2.
h1 = 0.1;
h2 = 0.1;

%% 3. Calculate reflection and transmission coefficients.
fluid.v = 1500;
fluid.density = 1000;
solid.v = 5900;
solid.vShear = 3150;
solid.density = 7850;
d = 12.4e-3;
model = MultiLayerModel(fluid, solid, fluid, d);

f = 0:0.5e3:600e3;
lambda = solid.v./f;

theta = 0.01:0.01:0.1;
%theta = 0;
%theta = 0.08;

xl = 0.52;
xh = 0.64;

nn = 2.^(10:15);

for n = nn
    x = linspace(0, 0.8, n);
    f = 0.5*x*solid.v/d;
    [V, W, k_hor, k_vert_L, debug] = fluidSolidFluid(f, theta, model);
end
legend(cellfun(@(x) num2str(x), num2cell(nn), 'UniformOutput', 0))
xlabel('\theta')
ylabel('Energy')



%% 5. Calculate the integral, phase shifting the wave and
%    multiplying the reflection/transmission coefficient
%% 6. Convolve with the excitation signal.