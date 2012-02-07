%% Plane wave normal incidence frequency response


%% Fluid-fluid-fluid model
fluid1 = struct('v', 1500, 'density', 1000);
fluid3 = fluid1;
layer = struct('v', 5900, 'density', 7850);
f = 100e3:1e3:1000e3;
%theta = 0.01:0.001:0.2;
theta = 0;
d = 12.4e-3;
[V, W] = fluidLayerReflectionCoefficient(f, theta, fluid3, layer, fluid1, d);
figure
subplot(211)
plot(f, real(V(:, 1)))
hold all
plot(f, real(W(:, 1)))
title('Real part')
legend('Reflection', 'Transmission')
subplot(212)
plot(f, imag(V(:, 1)))
hold all
plot(f, imag(W(:, 1)))
title('Imaginary part')