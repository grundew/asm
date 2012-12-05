f = 2000e3;
ntheta = 2.^20;
qmax = 1;

% Observation point
z = 10e-2;
x = 0;

%% Transducer specs
a = 9e-3;

%% Material parameters
rho_fluid = 1.5;
v_fluid = 350;
v_layer = 5850;
fluid1 = struct('v', v_fluid, 'density', rho_fluid);
fluid3 = fluid1;
layer = struct('v', v_layer, 'density', 7850, 'vShear', 3218);
d = 10.15e-3;

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, d);

tic
pquad = quadgk(@(xx) integrand(xx, f, a, x, z, model), -qmax, qmax,...
    'abstol', 1e-15, 'reltol', 0);
tquad = toc

tic
q = linspace(-qmax, qmax, ntheta);
dq = q(2) - q(1);
I = integrand(q, f, a, x, z, model);
ptrap = dq*(0.5*(I(1) + I(end)) + sum(I(2:end-1)));
ttrap = toc

pquad-ptrap