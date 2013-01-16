%% Test Levesque and Piche
freq = 200e3;
theta = linspace(1e-3, pi/2, 1e4);

%% Material parameters
% Water
rho_fluid = 1000;
c_F = 1500;
damping = 0;
% Air
% rho_fluid = 1;
% c_F = 350;
% v_layer = 5850 - i*1.25;
% v_shear = 3218 - i*1.25;
c_L = 5850;
c_S = 3218;
vShear_fluid = 0;

fluid1 = struct('v', c_F, 'density', rho_fluid, 'vShear', vShear_fluid);
fluid3 = fluid1;
layer = struct('v', c_L, 'density', 7850, 'vShear', c_S);
d = 10e-3;

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, d);

%% Run the model
[Ra, Ta] = analyticRTFast(freq, theta, model);
[Rb, Tb] = fluidSolidFluidReflectionCoefficient(freq, theta, model);

%% Plot the stuff
figure
subplot(211)
plot(theta, abs(Ra), '.', theta, abs(Rb), 'o')
title('Reflection coeffecient')
subplot(212)
plot(theta, unwrap(angle((Ra))), '.', theta, unwrap(angle((Rb))), 'o')
figure
subplot(211)
plot(theta, abs(Ta), '.', theta, abs(Tb), 'o')
title('Transmission coeffecient')
subplot(212)
plot(theta, unwrap(angle((Ta))), '.', theta, unwrap(angle((Tb))), 'o')