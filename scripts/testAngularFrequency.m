function testAngularFrequency()
% Test angular frequency

nfr = 5;
dfr = logspace(-1, -6, nfr);
f = 1000e3;

model = createmodel();

rtrapz = zeros(1, nfr);
figure
hold all
for i = 1:nfr
    % Calculate reflection coefficient for given angular resolution
    n_z = 1:-dfr(i):0;
    theta_z = acos(n_z);
    R = analyticRTFast(f, theta_z, model);
    rtrapz(i) = trapz(R, n_z);
    
    % Plot the absolute value of reflection coeffecient
    plot(theta_z, abs(R), '-')
end
legend(arrayfun(@(x) sprintf('%.2e', x), dfr, 'uni', 0))

% Plot the integrated reflection coefficient
figure
subplot(211)
plot(dfr, real(rtrapz))
subplot(212)
plot(dfr, imag(rtrapz))
end

function model = createmodel()
rho_fluid = 150;
v_fluid = 450;

v_layer = 5850;
fluid1 = struct('v', v_fluid, 'density', rho_fluid);
fluid3 = fluid1;
layer = struct('v', v_layer, 'density', 7850, 'vShear', 3162);
d = 25e-3;

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, d);
end