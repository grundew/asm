%% Plane wave normal incidence frequency response
fs = 2e6;
nfft = 2^27;
f = (0:nfft-1)*fs/nfft;

% Excitation pulse
tend = 50e-6;
t = (0:1/fs:tend)';

% Start and stop frequencies
f0 = 200e3;
f1 = 800e3;

% Window
wndw = gausswin(length(t));

% Real chirp
y = wndw.*chirp(t, f0, tend, f1, 'linear', 270);

% Pad with zeros
y = [zeros(100, 1); y; zeros(100, 1)];
tPulse = (0:length(y)-1)/fs';

% Make the waveform Analytic
xPulse = conj(hilbert(y));

%% Steel
v_solid = 5950;
v_solid_shear = 3218;

%% Gas steel water
% % A sort of natural gas (see aga10_testscript.m)
% rho_fluid1 = 110;
% v_fluid1 = 430;
% % Water
% rho_fluid2 = 1000;
% v_fluid2 = 1500;

%% Gas steel gas
% 100% Nitrogen (as in Kårstø)
rho_fluid1 = 156.01;
v_fluid1 = 408;

rho_fluid2 = rho_fluid1;
v_fluid2 = v_fluid1;

%% Fluid-fluid-fluid model
fluid1 = struct('v', v_fluid1, 'density', rho_fluid1);
fluid2 = struct('v', v_fluid2, 'density', rho_fluid2);
layer = struct('v', v_solid, 'density', 7850, 'vShear', v_solid_shear);

% Calculate angles
theta = 0;
theta_solid = asin(layer.v/fluid1.v*sin(theta));
theta_fluid2 = asin(fluid2.v/layer.v*sin(theta_solid));
theta_solid_shear = asin(layer.vShear/fluid1.v*sin(theta));
d = 30e-3;
fres = 0.5*5850/d;
dist = 0;

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid2, d);
[xR, xT, t] = planeWaveTimeSignal(model, xPulse, tPulse, theta, dist);

%% Estimate the decay rate per time (longitudenal)
Z = layer.density*layer.v/cos(theta_solid);
Z1 = fluid1.density*fluid1.v/cos(theta);
Z2 = fluid2.density*fluid2.v/cos(theta_fluid2);
R1 = (Z1 - Z)/(Z1 + Z);
R2 = (Z2 - Z)/(Z2 + Z);
alpha = -0.5*layer.v*log(R1*R2)/d;

yDecay = exp(-alpha*t);

%% Estimate the decay rate per time (shear)
Zshear = layer.density*layer.vShear/cos(theta_solid_shear);
R1shear = (Z1 - Zshear)/(Z1 + Zshear);
R2shear = (Z2 - Zshear)/(Z2 + Zshear);
alpha_shear = -0.5*layer.vShear*log(R1shear*R2shear)/d;

yDecay_shear = exp(-alpha_shear*t);


%% Plot the whole shebang
figure
plot(t, real(xR))
figure
a = 0.017;
semilogy(t*1e6, abs(xT), t*1e6, a*yDecay)
