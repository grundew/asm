clear;

%% Constants and parameters
r = 12e-3; % 12 cm
nfft = 2^12;
c = 1500;
rho = 1000;
f = 1500e3;
w = 2*pi*f;
k = w/c;

% Distance to reflector
h = 10e-2;

% q is sine theta
q = linspace(-0.7, 0.7, nfft);
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

%% Propagate the wave to all points on the x, z grid
tic();
p_i = propagateWave(X, Z, w, kx, kz, rho, V);
p_i = reshape(p_i, size(Z));
t1 = toc();
fprintf('Trapezoidal: %.3fs\n', t1);

%% Reflection coefficient
fluid.v = c;
fluid.density = rho;
solid.v = 5900;
solid.vShear = 3150;
solid.density = 7850;
d = 12.4e-3;
model = MultiLayerModel(fluid, solid, fluid, d);
[R, T] = fluidSolidFluidReflectionCoefficient(f, asin(q), model);

%% Propagate the reflected wave back
zr = fliplr(h:dz:2*h-dz);
[X, Zr] = meshgrid(x, zr);
tic();
p_r = propagateReflectedWave(X, Zr, w, kx, kz, rho, V, R);
p_r = reshape(p_r, size(Zr));
t1 = toc();
fprintf('Trapezoidal: %.3fs\n', t1);

%% Propagate the transmitted wave
zt = h:dz:2*h+d;
[Xt, Zt] = meshgrid(x, zt);
tic();
p_t = propagateReflectedWave(Xt, Zt, w, kx, kz, rho, V, T);
p_t = reshape(p_t, size(Zt));
t1 = toc();
fprintf('Trapezoidal: %.3fs\n', t1);

%% Plot the incoming pressure
figure
xx = [-fliplr(x), x];
pi = [fliplr(p_i), p_i];
imagesc(xx, z, abs(pi))
title('Incoming pressure')
clim = caxis();

%% Plot reflected pressure
figure
pr = [fliplr(p_r) p_r];
imagesc(xx, z, abs(pr))
title('Reflected pressure')
caxis(clim);

%% Plot total reflected pressure
figure
imagesc(xx, z, abs(pr+pi))
title('Total pressure')

%% Plot the complete pressure
figure
hold all
pt = [fliplr(p_t) p_t];

[XX, ZZ] = meshgrid([-fliplr(x) x], z);
[XXt, ZZt] = meshgrid([-fliplr(x) x], zt+d);
mesh(XX, ZZ, zeros(size(ZZ)), 10*log10(abs(pi+pr)))
mesh(XXt, ZZt, zeros(size(ZZt)), 10*log10(abs(pt)))

%% Plot the reflected and transmitted pressure
figure
hold all
title('Reflected and transmitted pressure')
pt = [fliplr(p_t) p_t];

[XX, ZZ] = meshgrid([-fliplr(x) x], z);
[XXt, ZZt] = meshgrid([-fliplr(x) x], zt+d);
mesh(XX, ZZ, zeros(size(ZZ)), 10*log10(abs(pr)))
mesh(XXt, ZZt, zeros(size(ZZt)), 10*log10(abs(pt)))