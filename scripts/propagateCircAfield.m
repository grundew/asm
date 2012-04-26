clear;

%% Constants and parameters
r = 12e-3; % 12 cm
nfft = 2^12;
c = 1500;
rho = 1000;
%f = [1000e3, 1500e3];
f = 100e3;
w = 2*pi*f;

% Distance to reflector
h = 10e-2;

% q is sine theta
q = linspace(-0.99, 0.99, nfft)';
qq = repmat(q, 1, length(f));
k = repmat(w/c, length(q), 1);
kx = k.*qq;
kz = k.*sqrt(1-qq.^2);

%% Circular aperture
V = angularPlaneWaveSpectrumPiston(r, c, q, f);
V = V/max(V);

%% Discretization of x, z
dx = 0.001;
dz = 0.001;
x = 0.001:dx:0.4;
z = 0.001:dz:h;
[X, Z] = meshgrid(x, z);

%% Propagate the wave to all points on the x, z grid
tic();
p_i = propagateWave(X, Z, f, q, kx, kz, rho, V);
t1 = toc();
fprintf('Trapezoidal: %.3fs\n', t1);

%% Plot the incoming pressure
xx = [-fliplr(x), x];
for i = 1:size(p_i, 3)
    figure
    pi = [fliplr(p_i(:, :, i)), p_i(:, :, i)];
    imagesc(xx, z, abs(pi))
    clim = caxis();
    title(sprintf('Incoming pressure. f = %d', f(i)))
    axis xy
end

%% Reflection coefficient
fluid.v = c;
fluid.density = rho;
solid.v = 5900;
solid.vShear = 3150;
solid.density = 7850;
d = 12.4e-3;
model = MultiLayerModel(fluid, solid, fluid, d);
tic();
[R, T] = fluidSolidFluidReflectionCoefficient(f, asin(q), model);
t1 = toc();
fprintf('Reflection coefficient: %.3fs\n', t1);

%% Propagate the reflected wave back
zr = fliplr(h:dz:2*h-dz);
[X, Zr] = meshgrid(x, zr);
tic();
p_r = propagateReflectedWave(X, Zr, f, q, kx, kz, rho, V, R);
t1 = toc();
fprintf('Trapezoidal: %.3fs\n', t1);

%% Plot reflected pressure
for i = 1:size(p_r, 3)
    pr = [fliplr(p_r(:, :, i)) p_r(:, :, i)];
    figure;
    imagesc(xx, z, db(abs(pr)));
    caxis(clim);
    title(sprintf('Reflected pressure. f = %d', f(i)));
end

%% Propagate the transmitted wave
zt = h:dz:2*h+d;
[Xt, Zt] = meshgrid(x, zt);
tic();
p_t = propagateReflectedWave(Xt, Zt, f, q, kx, kz, rho, V, T);
t1 = toc();
fprintf('Trapezoidal: %.3fs\n', t1);

%% Plot total reflected pressure
figure
imagesc(xx, z, db(abs(pr+pi)))
title('Total pressure')

%% Plot the complete pressure
figure
hold all
pt = [fliplr(p_t) p_t];

[XX, ZZ] = meshgrid([-fliplr(x) x], z);
[XXt, ZZt] = meshgrid([-fliplr(x) x], zt+d);
mesh(XX, ZZ, zeros(size(ZZ)), db(abs(pi+pr)))
mesh(XXt, ZZt, zeros(size(ZZt)), db(abs(pt)))
title('Total pressure and transmitted pressure')

%% Plot the reflected and transmitted pressure
figure
hold all
pt = [fliplr(p_t) p_t];

[XX, ZZ] = meshgrid([-fliplr(x) x], z);
[XXt, ZZt] = meshgrid([-fliplr(x) x], zt+d);
mesh(XX, ZZ, zeros(size(ZZ)), 10*log10(abs(pr)))
mesh(XXt, ZZt, zeros(size(ZZt)), 10*log10(abs(pt)))
title('Reflected and transmitted pressure')