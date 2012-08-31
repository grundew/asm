%% Fluid-fluid-fluid model
rho_fluid = 1000;
v_fluid = 1500;
v_layer = 5850;
v_shear = 3218;
fluid1 = struct('v', v_fluid, 'density', rho_fluid);
fluid3 = fluid1;
layer = struct('v', v_layer, 'density', 7850, 'vShear', v_shear);
d = 10.15e-3;
fres = 0.5*v_layer/d;

%% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, d);

% %% Matrix method. Single angle.
% theta = 0.01;
% f = 1e6:0.005e3:5e6;
% [V, W] = fluidSolidFluidReflectionCoefficient(f, theta, model);
% plotCoeffs(f/fres, theta, d, V, 'Reflection', 'Angle', '.-')
% plotCoeffs(f/fres, theta, d, W, 'Transmission', 'Angle', '.-')

%% Matrix method. Single frequency.

% At 2MHz there are problems
ntheta = 2^20;
theta = linspace(0, 0.5, ntheta);
f = 2000e3;
thresh = 1e-8;
[V, W] = fluidSolidFluidReflectionCoefficient(f, theta, model, thresh);

% Calculate all the stuff
% Horizontal wave number
klz = 2*pi*f/v_layer*sqrt((1 - v_layer/v_fluid*sin(theta)));
ksz = 2*pi*f/v_shear*sqrt((1 - v_shear/v_fluid*sin(theta)));
% Damping factor
alpha_L = imag(klz);
alpha_S = imag(ksz);
% Diagonal elements of the propagation vector
Cl = cos(klz*d);
Cs = cos(ksz*d);

plotCoeffs(theta, f/fres, d, V, 'Reflection', 'Angle', '.-')
plotCoeffs(theta, f/fres, d, W, 'Transmission', 'Angle', '.-')
plotCoeffs(theta, f/fres, d, klz, 'Horizontal wave number', 'Angle', '.-')
plotCoeffs(theta, f/fres, d, ksz, 'Horizontal wave number shear', 'Angle', '.-')

%% Plot the damping coefficient with the threshold
figure
semilogy(theta, exp(-alpha_L*d));
hold all
xl = xlim();
semilogy(theta, exp(-alpha_S*d));
semilogy(xl, thresh*ones(2, 1), 'r--')

%% Plot the diagonals in the propagation vector
figure
hold all
plot(theta, log10(abs(Cl)))
plot(theta, log10(abs(Cs)))

%% Investigate one point where it goes to hell
% The resulting value for W is about 65536
theta = 0.1731119;
f = 2000e3;
thresh = 1e-10;
[V, W, debug] = fluidSolidFluidReflectionCoefficient(f, theta, model, thresh);