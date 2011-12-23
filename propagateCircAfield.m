clear;

%% Constants and parameters
r = 12e-3; % 12 cm
nfft = 2^12;
c = 1500;
rho = 1000;
f = 100e3;
w = 2*pi*f;
k = w/c;

% Distance to reflector
h = 10e-2;

% q is sine theta
q = linspace(-0.999, 0.999, nfft);
kx = k*q;
kz = k*sqrt(1-q.^2);

%% Circular aperture
V = angularPlaneWaveSpectrumPiston(r, kx);

%% Discretization of x, z
dx = 0.001;
dz = 0.001;
x = 0.001:dx:0.4;
z = 0.001:dz:h;
[X, Z] = meshgrid(x, z);

%% Propagate the wave to all points on the x,z grid
tic();
p_i = propagateWave(Z, X, w, kx, kz, rho, V);
p_i = reshape(p_i, size(Z));
t1 = toc();
fprintf('Trapezoidal: %.3fs\n', t1);

%% Plotting
figure
xx = [-fliplr(x) x];
pi = abs([fliplr(p_i) p_i]);
imagesc(xx, z, abs(pi))
title('Incoming pressure')
clim = caxis();

%% Reflection coefficient
fluid.v = c;
fluid.density = rho;
solid.v = 5900;
solid.vShear = 3150;
solid.density = 7850;
d = 12.4e-3;
model = MultiLayerModel(fluid, solid, fluid, d);
[R, T] = fluidSolidFluid(f, asin(q), model);

%% Propagate the reflected wave back
zb = h+dz:dz:2*h;
[X, Zb] = meshgrid(x, zb);
tic();
p_r = propagateReflectedWave(Zb, X, w, kx, kz, rho, V, R);
p_r = reshape(p_r, size(Zb));
t1 = toc();
fprintf('Trapezoidal: %.3fs\n', t1);

%% Propagate the transmitted wave
zt = h+d:dz:2*h+d;
[Xt, Zt] = meshgrid(x, zt);
tic();
p_t = propagateReflectedWave(Zt, Xt, w, kx, kz, rho, V, T);
p_t = reshape(p_t, size(Zt));
t1 = toc();
fprintf('Trapezoidal: %.3fs\n', t1);

%% Plot reflected pressure
figure
pr = rot90([fliplr(p_r) p_r], 2);
imagesc(xx, z, abs(pr))
title('Reflected pressure')
caxis(clim);

%% Plot total reflected pressure
figure
imagesc(xx, z, abs(pr+pi))
title('Total pressure')

%% Plot the complete wave pressure
figure
hold all
pt = rot90([fliplr(p_t) p_t], 2);

[XX, ZZ] = meshgrid([-fliplr(x) x], z);
[XXt, ZZt] = meshgrid([-fliplr(x) x], zt);
mesh(XX, ZZ, zeros(size(ZZ)), 10*log10(abs(pi+pr)))
mesh(XXt, ZZt, zeros(size(ZZt)), 10*log10(abs(pt)))

%% Plot Various stuff
% figure
% hold all
% plot(q, abs(V))
% plot(q, abs(R))
% plot(q, abs(p_i(end, :)))
% legend('Spectrum of outgoing wave', 'Reflection coefficient', 'Pressure at reflection surface')