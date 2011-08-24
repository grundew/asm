clear all
%% Define parameters
theta = 10*pi/180; % Radians
ff = 100e3:0.01e3:1e6;
n = length(ff);
V = zeros(1, n);
W = zeros(1, n);

% Material props
solid = materials.MaterialFactory.produce('stainless steel');

c = 1500; % Speed of sound in water
rho_water = 1000;
rho_steel = 7850;
d = 12.4e-3;

% For keeping the determinant of A
deta = zeros(1, n);

for i = 1:n
    f = ff(i);
    w = 2*pi*f;
    
    % Wavenumbers
    
    % Length of wavenumber vector in fluid
    k = w/c;
    
    % Horizontal wavenumber (equal for all layers)
    K = k*sin(theta);
    
    % Vertical wavenumber in water
    k_z = k*cos(theta);
    
    % Length of wavenumber vector in the steel (S = shear, L = longitudenal)
    k_S = w/solid.vShear;
    k_L = w/solid.v;
    %k_S = sqrt(rho_steel*w^2/mu);
    %k_L = sqrt(rho_steel*w^2/(lame + 2*mu));
    
    % Horizontal part of wavenumber in steel
    k_z_S = sqrt(k_S^2 - K^2);
    k_z_L = sqrt(k_L^2 - K^2);
    
    %% Step 1:
    % Calculate input matrix
    input = inputMatrix(rho_water, w, k_z);
    
    %% Step 2:
    % Calculate the matrices related to each layer
    a = layerMatrix(rho_steel, w, k_z_S, k_z_L, K, k_S, d);
    
    %% Step 3:
    % Calculate the output matrix
    output = outputMatrix(rho_water, w, k_z);
    
    %% Step 4:
    % Calculate the products of all the similar contiguous matrices.
    
    % A = a(end, :, :);
    % for i = length(a)-1:-1:1
    %     A = a(i+1, :, :)*a(i, :, :);
    % end
    
    % Single layer
    A = a;
    deta(i) = det(A);
    
    %% Step 5:
    % Reduceing of the liquid-solid-liquid sandwiches.
    %BB = bDirectly(rho_steel, w, k_z_S, k_z_L, K, k_S, d);
    B = transformSolidMatrix(A);
    
    %% Calculate V and W
    G = output*B*input;
    
    V(i) = -G(2, 1)/G(2, 2);
    W(i) = G(1, 1) - G(1, 2)*G(2, 1)/G(2, 2);
end


%% Use Brekhovkikh to compare with
fluid = struct();                       % Define materials
layer = struct();                       % Define materials
fluid(1).cp = 1500;                     % Coupling medium, VOS
fluid(1).rho = 1000;                    % Coupling medium
cLayer = 1;

layer(cLayer).cp = solid.v;                % Solid layer, VOS
layer(cLayer).rho = solid.density;         % Solid layer
layer(cLayer).poisson = solid.poisson;     % Solid layer
layer(cLayer).thickness = d;               % Solid layer
fluid(2) = fluid(1);

msh = MultiShearModelPlateClass(ff, theta, layer, fluid);
msh.doAll();
R = msh.V;
T = msh.W;


%% Compare
figure,plot(ff, abs(V), '.', ff, abs(R), 'o')
figure,plot(ff, abs(W), '.', ff, abs(T), 'o')

figure
subplot(211)
plot(ff, real(V), '.', ff, real(R), 'o')
legend('Cervanka', 'Brekhovskikh')
title('Real part of Reflection coeff')
subplot(212)
plot(ff, imag(V), '.', ff, imag(R), 'o')
title('Imaginary part of Reflection coeff')

figure
plot(ff, abs(deta-1))
titlestr = sprintf('Determinant of a deveation from 1. Max %d', max(abs(deta-1)));
title(titlestr)