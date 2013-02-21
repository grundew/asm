function plotIntegrand_hei()

f = 200e3;
theta = linspace(0, pi/2, 2^12);
alphaLambda = 0;

% Distance from transducers to target
d1 = 10e-2;
d3 = 10e-2;

% Transducer radii
aRx = 6e-3;
aTx = 12e-3;

%% Material parameters
% rho_fluid = 1.5;
% v_fluid = 342.21;
rho_fluid = 1000;
v_fluid = 1500;

v_layer = 5850;
fluid1 = struct('v', v_fluid, 'density', rho_fluid);
fluid3 = fluid1;
layer = struct('v', v_layer, 'density', 7850, 'vShear', 3218);
d = 10e-3;

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, d);
deltax = 5e-2;

% I = orofinoIntegrand(theta, f, aRx, aTx,...
%     v_fluid, rho_fluid, d1, d3, model, alphaLambda);
I = shiftedTransducer(theta, f, aRx, aTx,...
   v_fluid, rho_fluid, d1, d3, model, alphaLambda, deltax);

subplot(211)
hold all
plot(180/pi*theta, abs(I)/max(abs(I)))
subplot(212)
hold all
plot(180/pi*theta, unwrap(angle(I)))