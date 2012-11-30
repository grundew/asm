%% Test one single frequency
debug = false;

%% Samplings stuff
ntheta = 2.^(10:2:18);
qmax = 1;

f0 = 500e3;
f1 = 1000e3;
nf = 10;
f = linspace(f0, f1, nf);

% Observation point
z = 10e-2;
x = 0;

%% Transducer specs
a = 9e-3;

%% Material parameters
rho_fluid = 1;
v_fluid = 350;
% v_layer = 5850 - i*1.25;
% v_shear = 3218 - i*1.25;
v_layer = 5850;
v_shear = 3218;
fluid1 = struct('v', v_fluid, 'density', rho_fluid);
fluid3 = fluid1;
layer = struct('v', v_layer, 'density', 7850, 'vShear', v_shear);
d = 10.15e-3;

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, d);
thresh = 1e-6;

p = zeros(nf, length(ntheta));

for j = 1:nf
    freq = f(j);
    %% Integrate over all angles for the point on the axis
    for i = 1:length(ntheta)
        % Uniform in q
        n = ntheta(i);
        %         if n == max(ntheta)
        %             debug = true;
        %         else
        %             debug = false;
        %         end
        
        % p(i) = integratePTrapezoidalStepWise(n, freq, x, z, model, a, qmax, nstep, false);
        % p(j, i) = integratePTrapezoidal(n, freq, x, z, model, a, qmax, debug);
        % p(j, i) = integratePHankelTransform(n, freq, x, z, model, a, qmax, debug);
        p(j, i) = integratePHankelTransformAdaptive(n, freq, x, z, model, a, qmax, debug);
        % p(j, i) = integratePHankelTransformSimpson(n, freq, x, z, model, a, qmax, debug);
    end
    
    %     %% Plot the pressure response for frequency freq
    %     figure(1)
    %     subplot(211)
    %     hold all
    %     loglog(ntheta, abs(real(p)), '.-')
    %     %loglog(ntheta, abs(real(p) - real(p(end))), '.-')
    %     ylabel('Real p')
    %     set(gca, 'YScale', 'log')
    %     set(gca, 'XScale', 'log')
    %     subplot(212)
    %     hold all
    %     loglog(ntheta, abs(imag(p)), '.-')
    %     ylabel('Imag p')
    %     xlabel('N_{\theta}')
    %     set(gca, 'YScale', 'log')
    %     set(gca, 'XScale', 'log')
end

%% Plot
figure
pfacit = repmat(p(:, end), 1, size(p, 2));
imagesc(ntheta, f, abs(p-pfacit)./abs(pfacit))
axis xy
colorbar
set(gca, 'Clim', [0 10])
xlabel('N_{\theta}')
ylabel('Frequency')