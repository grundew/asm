%% Test one single frequency

%% Samplings stuff
ntheta = 2.^(2:17);
thetamax = pi/2-0.01;

% Observation point
z = 10e-2;
x = 0;

%% Transducer specs
a = 12e-3;

%% Material parameters
rho_fluid = 1000;
v_fluid = 1500;
v_layer = 5850;
fluid1 = struct('v', v_fluid, 'density', rho_fluid);
fluid3 = fluid1;
layer = struct('v', v_layer, 'density', 7850, 'vShear', 3218);
d = 10.15e-3;

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, d);
thresh = 1e-9;

%% Integrate over all angles for the point on the axis
freq = 100e3;
p = zeros(length(ntheta), 1);

figure
subplot(211)
hold all
subplot(212)
hold all
for i = 1:length(ntheta)
    theta = linspace(-thetamax, thetamax, ntheta(i))';
    q = sin(theta);

    % Vertical and horizontal wave number
    kz = 2*pi*freq*cos(theta)/v_fluid;
    kx = 2*pi*freq*q/v_fluid;
    
    % Compute the factors in the integrand
    Phi = angularPlaneWaveSpectrumPiston(a, v_fluid, q, freq);
    [~, T] = fluidSolidFluidReflectionCoefficient(freq, theta, model, thresh);
    p(i) = propagateReflectedWave(x, z,...
        freq, q, kx, kz, rho_fluid, Phi, T);
    
    % Plot the integrands
    subplot(211)
    plot(theta, db(abs(T)/max(abs(T))))
    subplot(212)
    plot(theta, db(abs(Phi)/max(abs(Phi))), '.')
end

%% Plot the pressure response for frequency freq
figure
subplot(211)
plot(ntheta, real(p), '.-')
subplot(212)
plot(ntheta, imag(p), '.-')