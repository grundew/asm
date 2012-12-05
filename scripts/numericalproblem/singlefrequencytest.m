%% Test one single frequency
debug = false;
figure(1)

%% Samplings stuff
ntheta = 2.^(2:2:20);
qmax = 1;
freq = 2000e3;

% Observation point
z = 10e-2;
x = 1e-2;

%% Transducer specs
a = 9e-3;

%% Material parameters
% Water
rho_fluid = 1000;
c_F = 1500;
damping = 0;
% Air
% rho_fluid = 1;
% c_F = 350;
% v_layer = 5850 - i*1.25;
% v_shear = 3218 - i*1.25;
c_L = 5850;
c_S = 3218;
fluid1 = struct('v', c_F, 'density', rho_fluid);
fluid3 = fluid1;
layer = struct('v', c_L, 'density', 7850, 'vShear', c_S);
d = 10.15e-3;

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, d);
thresh = 1e-6;

%% Integrate over all angles for the point on the axis
p = zeros(length(ntheta), 1);
tid = zeros(length(ntheta), 1);
nstep = 200;

% % Calculate angles of refraction
% theta_L = asin(c_L/c_F*sin(theta)); 
% theta_S = asin(c_S/c_F*sin(theta));

for i = 1:length(ntheta)
    % Uniform in q
    n = ntheta(i);
    %     if n == max(ntheta)
    %         debug = true;
    %     else
    %         debug = false;
    %     end
    
    tic
    % p(i) = integratePTrapezoidalStepWise(n, freq, x, z, model, a, qmax, nstep, false);
    % p(i) = integratePTrapezoidal(n, freq, x, z, model, a, qmax, debug);
    p(i) = integratePHankelTransform(n, freq, x, z, model, a, qmax, debug);
    % p(i) = integratePHankelTransformSimpson(n, freq, x, z, model, a, qmax, debug);
    % p(i) = integratePHankelTransformAdaptive(n, freq, x, z, model, a, qmax, damping, debug);
    tid(i) = toc;
end

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

figure(2)
hold all
loglog(ntheta, tid)
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')