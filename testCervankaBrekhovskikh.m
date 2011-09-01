clear all
%% Define parameters
theta = 0.0159;
%theta = 1*pi/180; % Radians
ff = 100e3:0.5e3:1e6;

% Material props
solid = materials.MaterialFactory.produce('stainless steel');
fluidm = materials.MaterialFactory.produce('water');
d = 12.4e-3;
model = MultiLayerModel(fluidm, solid, fluidm, d);

[V, W] = fluidSolidFluid(ff, theta, model);

% Plotting

% % Reflection coefficient
% figure
% subplot(211)
% plot(ff, real(V))
% subplot(212)
% plot(ff, imag(V))
% 
% % Transmission coefficient
% figure
% subplot(211)
% plot(ff, real(W))
% subplot(212)
% plot(ff, imag(W))

%% Use Brekhovkikh to compare with
fluid = struct();                       % Define materials
layer = struct();                       % Define materials
fluid(1).cp = fluidm.v;                     % Coupling medium, VOS
fluid(1).rho = fluidm.density;                    % Coupling medium
cLayer = 1;

layer(cLayer).cp = solid.v;                % Solid layer, VOS
layer(cLayer).rho = solid.density;         % Solid layer
layer(cLayer).poisson = solid.poisson;     % Solid layer
layer(cLayer).thickness = d;               % Solid layer
fluid(2) = fluid(1);

msh = MultiShearModelPlateClass(ff, theta, layer, fluid);
msh.doAll();
R = msh.V;
T = msh.W;

%% Compare
figure,plot(ff, abs(V), '.', ff, abs(R), 'o')
figure,plot(ff, abs(W), '.', ff, abs(T), 'o')

figure
subplot(211)
plot(ff, real(V), '.', ff, real(R), 'o')
legend('Cervanka', 'Brekhovskikh')
title('Real part of Reflection coeff')
subplot(212)
plot(ff, imag(V), '.', ff, imag(R), 'o')
title('Imaginary part of Reflection coeff')

% figure
% plot(ff, abs(deta-1))
% titlestr = sprintf('Determinant of a deveation from 1. Max %d', max(abs(deta-1)));
% title(titlestr)