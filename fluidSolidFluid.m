function [V, W, k_hor, k_vert_L] = fluidSolidFluid(freq, theta_in, model, thresh)
% Modeling a water steel water system using method from Cervanka without
% taking into account over attenuated longitudenal waves.

if ~exist('thresh', 'var')
    % Threshold for where the algorithm stops propagating longitudenal
    % waves. Might be different for different thicknesses.
    thresh = 4.4e-15;
end

% Define parameters
nt = length(theta_in);
nf = length(freq);
V = zeros(nf, nt);
W = zeros(nf, nt);

k_hor = zeros(nf, nt);
k_vert_L = zeros(nf, nt);

% Speed of sounds
c_front = model.fluid(1).v;
c_back = model.fluid(2).v;
vShear = model.solid.vShear;
vLong = model.solid.v;

% Densities
rho_fluidFront = model.fluid(1).density;
rho_fluidBack = model.fluid(2).density;
rho_solid = model.solid.density;

% Thickness of plate
d = model.thickness;

% Absorbsion
% alpha_front = 1000;
% alpha_back = 1;
% alpha_L = 10;
% alpha_S = 10;


for i = 1:nf
    
    f = freq(i);
    w = 2*pi*f;
    
    % Wavenumbers
    
    % Length of wavenumber vector in the steel (S = shear, L = longitudenal)
    k_S = w/vShear;
    k_L = w/vLong;
        
    % Length of wavenumber vector in fluids
    k_front = w/c_front;
    k_back = w/c_back;
    for j = 1:nt
        theta = theta_in(j);
        % Horizontal wavenumber (equal for all layers)
        K = k_front*sin(theta);
        k_hor(i, j) = K;
        
        % Vertical wavenumber in front fluid
        k_z_front = k_front*cos(theta);
        
        % Horizontal part of wavenumber in steel
        k_z_S = sqrt(k_S^2 - K^2);
        k_z_L = sqrt(k_L^2 - K^2);
        k_vert_L(i, j) = k_z_L;
                
        % Step 1:
        % Calculate input matrix
        input = inputMatrix(rho_fluidFront, w, k_z_front);
        
        % Step 2:
        % Calculate the matrices related to each layer
        % a = layerMatrix(rho_steel, w, k_z_S, k_z_L, K, k_S, d);
        if theta == 0
            % If angle is zero the shear waves in solid is neglected.
            B = fluidLayerMatrix(rho_solid, w, k_z_L, d);
        else
            B = fluidSolidFluidLayer(rho_solid, w, k_z_S, k_z_L, K, k_S, d, thresh);
        end
        
        % Step 3:
        % Calculate the output matrix.
        k_z_back = sqrt(k_back^2 - K^2);
        output = outputMatrix(rho_fluidBack, w, k_z_back);
        
        % Calculate V and W
        G = output*B*input;
        
        V(i, j) = -G(2, 1)/G(2, 2);
        W(i, j) = G(1, 1) - G(1, 2)*G(2, 1)/G(2, 2);
    end
    
end