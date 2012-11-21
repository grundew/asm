%% Test one single frequency

%% Samplings stuff
ntheta = 2.^(11:19);
thetamax = 0.1;
freq = 2000e3;

% Observation point
z = 10e-2;
x = 0;

%% Transducer specs
a = 12e-3;

%% Material parameters
rho_fluid = 1;
v_fluid = 350;
v_layer = 5850;
fluid1 = struct('v', v_fluid, 'density', rho_fluid);
fluid3 = fluid1;
layer = struct('v', v_layer, 'density', 7850, 'vShear', 3218);
d = 10.15e-3;

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, d);
thresh = 1e-6;

%% Integrate over all angles for the point on the axis
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
    Phi = planePistonPressureAngularSpectrum(z, kx, kz, a, v_fluid, rho_fluid);
    [R, T] = analyticRTFast(freq, theta, model);

    Ht = Phi.*T.*exp(1i*kz*z);
    
    E = exp(1i*kx.*x);
    p = trapz(q, Ht.*E, 1);
    
    %     p(i) = propagateReflectedWave(x, z,...
    %         freq, q, kx, kz, rho_fluid, Phi, T);

    % Plot the integrands
    ax(1) = subplot(211);
    plot(theta, db(abs(R)))
    ax(2) = subplot(212);
    plot(theta, db(abs(Phi)/max(abs(Phi))), '.')
    linkaxes(ax, 'x')
end

%% Plot the pressure response for frequency freq
figure
subplot(211)
plot(ntheta, real(p), '.-')
subplot(212)
plot(ntheta, imag(p), '.-')