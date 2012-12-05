%% Test one single frequency

%% Samplings stuff
ntheta = 2^12;
qmax = 0.5;
freq = 2000e3;

% Observation point
z = 10e-2;
x = 0;

%% Transducer specs
a = 9e-3;

%% Angular spectrum dependenc on damping
% With damping
rho_fluid = 1000;
v_fluid = 1500;
v_layer = 5850 - 1i;
vshear = 3218 - 1i;
fluid1 = struct('v', v_fluid, 'density', rho_fluid);
fluid3 = fluid1;
layer = struct('v', v_layer, 'density', 7850, 'vShear', vshear);
d = 10e-3;

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, d);
q = linspace(-qmax, qmax, ntheta)';
[~, Tdamped] = analyticRTFast(freq, asin(q), model);

% Without damping
rho_fluid = 1000;
v_fluid = 1500;
v_layer = 5850;
vshear = 3218;
fluid1 = struct('v', v_fluid, 'density', rho_fluid);
fluid3 = fluid1;
layer = struct('v', v_layer, 'density', 7850, 'vShear', vshear);
d = 10e-3;

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, d);
q = linspace(-qmax, qmax, ntheta)';
[~, T] = analyticRTFast(freq, asin(q), model);

% Plot comparison
figure
plot(q, db(abs(T)), '.', q, db(abs(Tdamped)), 'o')
xlabel('sin \theta')
ylabel('abs(T)')
legend('No damping', 'Damping')

%% Frequency spectrum dependenc on damping at normal incidence
freq = 50e3:1e3:1e6;
theta = 0;

% With damping
rho_fluid = 1000;
v_fluid = 1500;
v_layer = 5850 - 1i;
vshear = 3218 - 1i;
fluid1 = struct('v', v_fluid, 'density', rho_fluid);
fluid3 = fluid1;
layer = struct('v', v_layer, 'density', 7850, 'vShear', vshear);
d = 10e-3;

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, d);
[~, Tdamped] = analyticRTFast(freq, theta, model);

% Without damping
rho_fluid = 1000;
v_fluid = 1500;
v_layer = 5850;
vshear = 3218;
fluid1 = struct('v', v_fluid, 'density', rho_fluid);
fluid3 = fluid1;
layer = struct('v', v_layer, 'density', 7850, 'vShear', vshear);
d = 10e-3;

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, d);
[~, T] = analyticRTFast(freq, theta, model);

% Plot comparison
figure
plot(freq, db(abs(T)), '.', freq, db(abs(Tdamped)), 'o')
xlabel('Frequency')
ylabel('abs(T)')
legend('No damping', 'Damping')