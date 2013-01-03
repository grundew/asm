%% Samplings stuff
nf = 2^12;
f = linspace(0, 1000e3, nf);
theta = 0.01;
q = sin(theta);
ntest = 3;

% Physical parameters
rho_fluid = 500;
cfluid = linspace(1000, 3000, ntest);
rho_plate = 7850;
c_L = 5850;
c_S = 3218;

figure(1)
hold all
for i = 1:ntest
    c_F = cfluid(i);
    fluid1 = struct('v', c_F, 'density', rho_fluid);
    fluid3 = fluid1;
    layer = struct('v', c_L, 'density', rho_plate, 'vShear', c_S);
    d = 10e-3;

    % Helper parameters
    theta_2 = asin(c_L/c_F*sin(theta));
    Z1 = c_F*rho_fluid/cos(theta);
    Z2 = c_L*rho_plate/cos(theta_2);
    fres = 0.5*c_L/d;

    % Scaling
    scalingx = 1/fres;
    % scalingy = Z2/Z1;
    scalingy = 1;

    % Fluid-solid-fluid model
    model = MultiLayerModel(fluid1, layer, fluid3, d);
    T = transmissionCoefficientAnalytical(f, q, model);
    
    %% Plot the stuff
    plot(f*scalingx, abs(T)*scalingy)
end
