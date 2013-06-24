function V = fluidSolidSolidReflectionCoefficient(freq, theta_in, model, thresh)
% Modeling a fluid-solid-solid system using method from Cervanka without
% taking into account over attenuated longitudenal waves.

if ~exist('thresh', 'var')
    % Threshold for where the algorithm stops propagating longitudenal
    % waves. Might be different for different thicknesses.
    thresh = 4.5e-15;
end

% Define parameters
nt = length(theta_in);
nf = length(freq);
V = zeros(nt, nf);

% Speed of sounds
c_front = model.fluid(1).v;
c_back = model.fluid(2).v;
cShear_back = model.fluid(2).vShear;
vShear = model.solid.vShear;
vLong = model.solid.v;

% Densities
rho_fluidFront = model.fluid(1).density;
rho_solidBack = model.fluid(2).density;
rho_solid = model.solid.density;

% Thickness of plate
d = model.thickness;

for i = 1:nf
    
    f = freq(i);
    w = 2*pi*f;
    
    % Length of wavenumber vector in the steel (S = shear, L = longitudenal)
    k_S = w/vShear;
    k_L = w/vLong;
    
    % Length of wavenumber vector in input fluid
    k_front = w/c_front;
    
    % Length of wavenumber vectors in output solid
    % k_L_back = w/c_back;
    k_S_back = w/cShear_back;
    
    for j = 1:nt
        theta = theta_in(j);
        % Horizontal wavenumber (equal for all layers)
        K = k_front*sin(theta);
        
        % Vertical wavenumber in front fluid
        k_z_front = k_front*cos(theta);
        
        % Horizontal part of wavenumber in solid layer
        k_z_S = sqrt(k_S^2 - K^2);
        k_z_L = sqrt(k_L^2 - K^2);

        % Step 1:
        % Calculate input matrix
        
        % Dimensions are wrong
        % 2 x 2
        rhow2m = -rho_fluidFront*w^2;
        input = [k_z_front, -k_z_front;...
            rhow2m, rhow2m];
        
        % Step 2:
        % Calculate the matrices related to each layer
        % a = layerMatrix(rho_steel, w, k_z_S, k_z_L, K, k_S, d);
        if theta == 0
            % If angle is zero the shear waves in solid is neglected.
            
            % What to do here?
            % B is 2 x 2
            B = fluidLayerMatrix(rho_solid, w, k_z_L, d);
        else
            % 4 x 4
            A = solidLayerMatrix(rho_solid, w, k_z_S, k_z_L, K, k_S, d);
        end
        
        % Step 3:
        % Calculate the output matrix.
        k_z_L_back = sqrt(k_L_back^2 - K^2);
        k_z_S_back = sqrt(k_S_back^2 - K^2);
        
        % 4 x 4
        output = outputSolidMatrix(rho_solidBack, w, k_z_S_back, k_z_L_back, K, k_S_back);
        
        % Calculate V
        H = output*A;
        G = input;
        
        
        V(j, i) = 1;
        V(j, i) = -G(2, 1)/G(2, 2);
    end
    
end

fzero = freq==0;
V(:, fzero) = zeros(nt, nnz(fzero));
end