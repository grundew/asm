%% Define parameters
theta = 0.1; % Radians
f = 100e3;
w = 2*pi*f;

% Material props
c = 1450; % Speed of sound in water
rho_water = 1000;
rho_steel = 7810;
d = 13.4e-3;
mu = 7.79e10; % Lame parameters
lame = 6.7625e10;

% Wavenumbers
% Horizontal wavenumber (equal for all layers)
K = w/c*sin(theta);

% Vertical wavenumber in water
k_z = w/c*cos(theta);

% Length of wavenumber in the steel
k_S = sqrt(rho_steel*w^2/mu);
k_L = sqrt(rho_steel*w^2/(lame + 2*mu));

% Horizontal part of wavenumber in steel
k_z_S = sqrt(k_S^2 - K^2);
k_z_L = sqrt(k_L^2 - K^2);

%% Step 1:
% Calculate input matrix
input = inputMatrix(rho_water, w, k_z);

%% Step 2:
% Calculate the matrices related to each layer
a = layerMatrix(rho_steel, w, k_z_S, k_z_L, K, d, mu);

%% Step 3:
% Calculate the output matrix
output = outputMatrix(rho_water, w, k_z_L);

%% Step 4:
% Calculate the products of all the similar contiguous matrices.

% A = a(end, :, :);
% for i = length(a)-1:-1:1
%     A = a(i+1, :, :)*a(i, :, :);
% end

% Single layer
A = a;

%% Step 5:
% Reduceing of the liquid-solid-liquid sandwiches.
B = transformSolidMatrix(A);

%% Calculate V and W
G = output*B*input;

V = -G(2, 1);
W = G(1, 1) - G(1, 2)*G(2, 1);