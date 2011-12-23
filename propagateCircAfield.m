clear;

%% Constants and parameters
r = 12e-3; % 12 cm
nfft = 2^12;
c = 1500;
rho = 1000;
f = 100e3;
w = 2*pi*f;
k = w/c;
% q is sine theta
q = linspace(-0.999, 0.999, nfft);
kx = k*q;
kz = k*sqrt(1-q.^2);

%% Circular aperture
V = angularPlaneWaveSpectrumPiston(r, kx);

%% Discretization of x, z
z = 0.01:0.001:0.5;
x = 0.001:0.001:0.4;

%% Propagate the wave to all points on the x,z grid
tic();
p_i = propagateWave(z, x, w, kx, kz, rho, V);
t1 = toc();

fprintf('Trapezoidal: %.3fs\n', t1);

%% Reflect and propgate back
% Ideal reflector
R = ones(size(q));

%% Plotting
figure
title('Incoming pressure')
imagesc([-fliplr(x) x], z, abs([fliplr(p_i) p_i]))

%% Propagate back
