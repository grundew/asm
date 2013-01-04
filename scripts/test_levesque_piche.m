%% Test Levesque and Piche
freq = 200e3;
theta = linspace(1e-3, pi/2, 10);

%% Material parameters
% Water
rho_fluid = 1.5;
c_F = 340;
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
[Rp, Tp, Rs, Ts] = reflectionCoefficientLayeredMedia(freq, theta, model);


%% Plot the stuff
plot(theta, abs(Rp))