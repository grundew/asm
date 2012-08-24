%% Fluid-fluid-fluid model
rho_fluid = 1000;
v_fluid = 1500;
v_layer = 5850;
fluid1 = struct('v', v_fluid, 'density', rho_fluid);
fluid3 = fluid1;
layer = struct('v', v_layer, 'density', 7850, 'vShear', 3218);
theta = 0;
d = 10.15e-3;

%% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, d);
nfft = 2^12;
fs = 2e6;
f = fftshift((-nfft/2:nfft/2-1)*fs/nfft);

%% Analytical expression
[V, W] = analyticRT(f, theta, model);

%% Matrix method
tic
[VV, WW] = fluidSolidFluidReflectionCoefficient(f, theta, model);
toc

%% Plot both the reflection coefficients
figure
subplot(211),plot(f, real(VV), '.', f, real(V), 'o')
subplot(212),plot(f, imag(VV), '.', f, imag(V), 'o')

%% Plot both the transmission coefficients
figure
subplot(211),plot(f, real(WW), '.', f, real(W), 'o')
subplot(212),plot(f, imag(WW), '.', f, imag(W), 'o')

%% Calculate for several angles
ntheta = 256;
theta = linspace(0, 0.9, ntheta);
f = 100e3:1e3:900e3;
profile on
[R, T] = fluidSolidFluidReflectionCoefficient(f, theta, model);
profile viewer

%% Plot the beast
figure
subplot(211)
imagesc(f, theta, real(R));
axis xy
subplot(212)
imagesc(f, theta, imag(R));
axis xy