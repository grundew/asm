a = 12e-3;
h = 10e-2;
fluid.v = 1500;
fluid.density = 1000;
solid.v = 5900;
solid.vShear = 3150;
solid.density = 7850;
d = 10.15e-3;
model = MultiLayerModel(fluid, solid, fluid, d);
nq = 2^10;

% %% Test the geometry
% f = 100e3;
% wni = WaveNumberIntegration2(a, h, f, model, nq);
% 
% %% Discretization of x, z
% dx = 0.001;
% dz = 0.001;
% x = 0.001:dx:0.4;
% z = 0:dz:h-dz;
% [X, Z] = meshgrid(x, z);
% tic
% P = calculateReflectedPressure(wni, X, Z);
% toc
% 
% %% Plot reflected pressure
% xx = [-fliplr(x), x];
% pr = [fliplr(P) P];
% figure
% imagesc(xx, z, db(abs(pr)))
% title('Reflected pressure')
% 
% %% Discretization of x, z for transmitted area
% zt = h:dz:2*h+d;
% [Xt, Zt] = meshgrid(x, zt);
% Pt = wni.calculateTransmittedPressure(Xt, Zt);
% 
% %% Plot transmitted pressure
% xx = [-fliplr(x), x];
% pt = [fliplr(Pt) Pt];
% figure
% imagesc(xx, zt, db(abs(pt)))
% title('Transmitted pressure')
% 
% %% Plot transmitted and reflected pressure
% [XX, ZZ] = meshgrid([-fliplr(x) x], z);
% [XXt, ZZt] = meshgrid([-fliplr(x) x], zt+d);
% mesh(XX, ZZ, zeros(size(ZZ)), db(abs(pr)))
% mesh(XXt, ZZt, zeros(size(ZZt)), db(abs(pt)))
% title('Reflected and transmitted pressure')

%% Test time signal
% Chirp pulse
fs = 50e6;
nfft = 2^14;
tend = 50.7e-6;
ty_tmp = (0:1/fs:tend)';
f0 = 50e3;
f1 = 600e3;

alpha = (f1-f0)/tend;
wndw = gausswin(length(ty_tmp));
y = wndw.*1j.*exp(1j*2*pi*alpha*ty_tmp.^2/2);
y = [zeros(100, 1); y; zeros(100, 1)];
ty = 1/fs*(0:length(y)-1)';
Y = ifft(y, nfft);

f = (0:nfft-1)*fs/nfft;
% Real chirp
% wndw = gausswin(length(t));
% y = wndw.*chirp(t, f0, tend, f1, 'linear', 270);
% y = [zeros(100, 1); y; zeros(100, 1)];
% t = (0:length(y)-1)/fs';
%y = chirp(t, f0, t(end), f1, 'linear', 270);

% Filter to remove the DC-component
% n = 200;
% wn = [10e3, 2000e3]./fs/2;
% b = fir1(n, wn);
% y = filtfilt(b, 1, y);

%% Plot pulse
figure
plot(ty, real(y))

%% Calculate time signal
wni = WaveNumberIntegration2(a, h, [], model, nq);

nx = 10;
% [x, tx] = calculateReflectedPressureRx(wni, a, h, nx, y, fs, f);
[x, tx] = calculateReflectedPressureRx(wni, a, h, nx, Y.', fs, f);

%% Plot the time signal
figure
plot(tx, x)
