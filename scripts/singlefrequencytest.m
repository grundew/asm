%% Test one single frequency

%% Samplings stuff
ntheta = 2.^(8:2:30);
qmax = 1;
freq = 500e3;

% Observation point
z = 10e-2;
x = 0;

%% Transducer specs
a = 9e-3;

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
pp = zeros(length(ntheta), 1);
nstep = 20;

tic
for i = 1:length(ntheta)
    % Uniform in q
    n = ntheta(i);
    % p(i) = integratePTrapezoidal(n, freq, x, z, model, a, qmax);
    p(i) = integratePTrapezoidalStepWise(n, freq, x, z, model, a, qmax, nstep);
end
toc

%% Plot the pressure response for frequency freq
figure(1)
subplot(211)
hold all
loglog(ntheta, abs(real(p)), '.-')
%loglog(ntheta, abs(real(p) - real(p(end))), '.-')
ylabel('Real p')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
subplot(212)
hold all
loglog(ntheta, abs(imag(p)), '.-')
ylabel('Imag p')
xlabel('N_{\theta}')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

% % Plot the integrand
% figure
% subplot(211)
% plot(q, db(abs(Ht)))
% subplot(212)
% plot(q, unwrap(angle(Ht)))