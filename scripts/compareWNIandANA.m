fs = 2e6;
nfft = 2^12;
f = fftshift((-nfft/2:nfft/2-1)*fs/nfft);
% f = -200e3;

fluid1 = struct('v', 1500, 'density', 1000);
% fluid1 = struct('v', 340, 'density', 1000);
fluid3 = fluid1;
layer = struct('v', 5850, 'density', 7850, 'vShear', 3218);
theta = 0.1;
%theta = 1e-4;
d = 10.15e-3;

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, d);

tic
[V, W, K, tL, tS, kL, kS, kzL, kzS] = analyticRT(f, theta, model);
toc
tic
[R, T, KK, ttL, ttS, kkL, kkS, kkzL, kkzS] = fluidSolidFluidReflectionCoefficient(f, theta, model);
toc

% Calculate horizontal and total wave vector components


%% Plot the scheisen
figure
plot(f, K, '.', f, KK, 'o')
title('Horizontal wave number')
figure
plot(f, tL, '.', f, ttL, 'o')
title('Theta longitudenal')
figure
plot(f, tS, '.', f, ttS, 'o')
title('Theta shear')
figure
plot(f, kL, '.', f, kkL, 'o')
title('k longitudenal')
figure
plot(f, kS, '.', f, kkS, 'o')
title('k shear')
figure
subplot(211)
plot(f, real(V), '.', f, real(R), 'o')
title('Reflection coefficient')
subplot(212)
plot(f, imag(V), '.', f, imag(R), 'o')
figure
subplot(211)
plot(f, real(W), '.', f, real(T), 'o')
title('Transmission coefficient')
subplot(212)
plot(f, imag(W), '.', f, imag(T), 'o')
figure
plot(f, sqrt(K.^2 + kzL.^2), '.', f, sqrt(KK.^2 + kkzL.^2), 'o')
hold all,plot(f, abs(2*pi*f/model.solid.v), 'rd')
title('Total wave number longitudenal')
figure
plot(f, kzL, '.', f, kkzL, 'o')
title('Vertical wavenumber longitudenal')
figure
plot(f, kzS, '.', f, kkzS, 'o')
title('Vertical wavenumber shear')
figure
plot(f, sqrt(K.^2 + kzS.^2), '.', f, sqrt(KK.^2 + kkzS.^2), 'o')
hold all,plot(f, abs(2*pi*f/model.solid.vShear), 'rd')
title('Total wave number shear')
