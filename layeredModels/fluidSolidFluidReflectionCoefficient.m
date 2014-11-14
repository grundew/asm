function V = fluidSolidFluidReflectionCoefficient(freq, theta_in, model, thresh)
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
W = zeros(nt, nf);

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

alpha_L = model.alphaLambda*log(10.^(abs(freq)/model.solid.v/20));
alpha_S = (alpha_L/2)*(vLong/vShear)^3;

% Equation 5 & 6 from Ref(2)
vLong = vLong./(1 + 1i*alpha_L*vLong./freq/2/pi);
vShear = vShear./(1 + 1i*alpha_S*vShear./freq/2/pi);

for i = 1:nf
    
    f = freq(i);
    w = 2*pi*f;
    
    % Wavenumbers
    
    % Length of wavenumber vector in the steel (S = shear, L = longitudenal)
    k_S = w./vShear(i);
    k_L = w./vLong(i);
    
    % Length of wavenumber vector in fluids
    k_front = w/c_front;
    k_back = w/c_back;
    for j = 1:nt
        theta = theta_in(j);
        % Horizontal wavenumber (equal for all layers)
        K = k_front*sin(theta);
        
        % Vertical wavenumber in front fluid
        k_z_front = k_front*cos(theta);
        
        % Horizontal part of wavenumber in steel
        k_z_S = sqrt(k_S.^2 - K^2);
        k_z_L = sqrt(k_L.^2 - K^2);

        % Step 1:
        % Calculate input matrix
        rhow2m = -rho_fluidFront*w^2;
        input = [k_z_front, -k_z_front;...
            rhow2m, rhow2m];
        
        % Step 2:
        % Calculate the matrices related to each layer
        % a = layerMatrix(rho_steel, w, k_z_S, k_z_L, K, k_S, d);
        if theta == 0
            % If angle is zero the shear waves in solid is neglected.
            B = fluidLayerMatrix(rho_solid, w, k_z_L, d);
        else
            B = fluidSolidFluidLayerMatrix(rho_solid, w, k_z_S, k_z_L, K, k_S, d, thresh);
        end
        
        % Step 3:
        % Calculate the output matrix.
        k_z_back = sqrt(k_back.^2 - K.*2);
        rhow2 = rho_fluidBack*w^2;
        output = [0, -1/rhow2;...
            -1, -k_z_back/rhow2];
        % Calculate V and W
        G = output*B*input;
        
        if isequal(G, zeros(2, 2))
            V(j, i) = 1;
            W(j, i) = 0;
        else
            V(j, i) = -G(2, 1)/G(2, 2);
            W(j, i) = G(1, 1) - G(1, 2)*G(2, 1)/G(2, 2);
        end

    end
    
end

fzero = freq==0;
V(:, fzero) = zeros(nt, nnz(fzero));
end