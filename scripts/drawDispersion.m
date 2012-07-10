%% Material parameters
fluid1 = struct('v', 1500, 'density', 1000);
fluid3 = fluid1;
layer = struct('v', 5850, 'density', 7850, 'vShear', 2000);
nt = 100;
theta = linspace(0, 1, nt);
%theta = 1e-4;
d = 10e-3;

nf = 2000;
fres = layer.v/2/d;
f =  linspace(0, 5*fres, nf);
% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, d);

%% Solid layer immersed in a fluid
%% Calculate V and W
V = analyticRT(f, theta, model);
%% Fluid layer immersed in a fluid
%% Calculate V and W
% Rf =  fluidLayerReflectionCoefficient(f, theta, fluid1, layer, fluid1, d);

%% Plot
figure
imagesc(f/fres, theta, abs(V)')
axis xy
xlabel('f/fres')
ylabel('\theta')
title(sprintf('Solid layer immersed in fluid. cL/cP = %f', layer.v/layer.vShear))
colormap(gray)

% %% Plot
% figure
% imagesc(f/fres, theta, abs(Rf)')
% axis xy
% xlabel('f/fres')
% ylabel('\theta')
% title('Fluid layer immersed in fluid')
% colormap(gray)