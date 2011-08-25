function [V, W, k_hor, k_vert_L, alpha_L] = waterSteelWaterC1(freq, theta_in, d, thresh)
% Modeling a water steel water system using method from Cervanka without
% taking into account over attenuated longitudenal waves.

% Define parameters
nt = length(theta_in);
nf = length(freq);
V = zeros(nf, nt);
W = zeros(nf, nt);

k_hor = zeros(nf, nt);
k_vert_L = zeros(nf, nt);

% Material props
solid.v = 5348.6;
solid.vShear = 3158.2;

c = 1500; % Speed of sound in water
rho_water = 1000;
rho_steel = 7850;

for i = 1:nf
    
    f = freq(i);
    w = 2*pi*f;
    
    % Wavenumbers
    
    % Length of wavenumber vector in the steel (S = shear, L = longitudenal)
    k_S = w/solid.vShear;
    k_L = w/solid.v;
        
    % Length of wavenumber vector in fluid
    k = w/c;
    
    for j = 1:nt
        theta = theta_in(j);
        % Horizontal wavenumber (equal for all layers)
        K = k*sin(theta);
        k_hor(i, j) = K;
        
        % Vertical wavenumber in water
        k_z = k*cos(theta);
        
        % Horizontal part of wavenumber in steel
        k_z_S = sqrt(k_S^2 - K^2);
        k_z_L = sqrt(k_L^2 - K^2);
        k_vert_L(i, j) = k_z_L;
                
        % Step 1:
        % Calculate input matrix
        input = inputMatrix(rho_water, w, k_z);
    
        % Step 2:
        % Calculate the matrices related to each layer
        %a = layerMatrix(rho_steel, w, k_z_S, k_z_L, K, k_S, d);
        if theta == 0
            % If angle is zero the shear waves in solid is neglected.
            B = fluidLayerMatrix(rho_steel, w, k_z_L);
        else
            B = fluidSolidFluidLayer(rho_steel, w, k_z_S, k_z_L, K, k_S, d, thresh);
        end

        % Step 3:
        % Calculate the output matrix
        output = outputMatrix(rho_water, w, k_z);
        
        % Calculate V and W
        G = output*B*input;
        
        V(i, j) = -G(2, 1)/G(2, 2);
        W(i, j) = G(1, 1) - G(1, 2)*G(2, 1)/G(2, 2);
    end
    
end