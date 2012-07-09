function [V, W, KK, thetal, thetas, KL, KS, k_vert_L, k_vert_S] = fluidSolidFluidReflectionCoefficient(freq, theta_in, model, thresh)
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

V = zeros(nt, nf);
W = zeros(nt, nf);

k_hor = zeros(nt, nf);
k_vert_L = zeros(nt, nf);
k_vert_S = zeros(nt, nf);

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

% Wave numbers to output
thetal = zeros(nf, nt);
thetas = zeros(nf, nt);

KL = zeros(nf, 1);
KS = zeros(nf, 1);

% Absorbsion
% alpha_front = 1000;
% alpha_back = 1;
% alpha_L = 10;
% alpha_S = 10;

for i = 1:nf
    
    f = freq(i);
    w = 2*pi*f;
    
    %     if f == 0
    %         V(:, i) = ones(nt, 1);
    %         W(:, i) = zeros(nt, 1);
    %         continue
    %     end
    % Wavenumbers
    
    % Length of wavenumber vector in the steel (S = shear, L = longitudenal)
    k_S = w/vShear;KS(i) = k_S;
    k_L = w/vLong;KL(i) = k_L;
        
    % Length of wavenumber vector in fluids
    k_front = w/c_front;
    k_back = w/c_back;
    for j = 1:nt
        theta = theta_in(j);
        % Horizontal wavenumber (equal for all layers)
        K = k_front*sin(theta);
        k_hor(j, i) = K;
        
        % Vertical wavenumber in front fluid
        k_z_front = k_front*cos(theta);
        
        % Horizontal part of wavenumber in steel
        k_z_S = sqrt(k_S^2 - K^2);
        k_z_L = sqrt(k_L^2 - K^2);
        % k_z_S = k_S*asin(K/k_L);
        % k_z_L = k_L*asin(K/k_L);

        k_vert_L(j, i) = k_z_L;
        k_vert_S(j, i) = k_z_S;
        % Debug variables
        thetal(i, j) = asin(K/k_L);
        thetas(i, j) = asin(K/k_S);
        
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
            B = fluidSolidFluidLayerMatrix(rho_solid, w, k_z_S, k_z_L, K, k_S, d, thresh);
        end
        
        % Step 3:
        % Calculate the output matrix.
        k_z_back = k_back*cos(theta_in);
        output = outputMatrix(rho_fluidBack, w, k_z_back);
        
        % Calculate V and W
        G = output*B*input;
        
        V(j, i) = -G(2, 1)/G(2, 2);
        W(j, i) = G(1, 1) - G(1, 2)*G(2, 1)/G(2, 2);
    end
    
end

KK = k_hor;

fzero = freq==0;
V(:, fzero) = zeros(nt, nnz(fzero));
W(:, fzero) = ones(nt, nnz(fzero)) + 1i*zeros(nt, nnz(fzero));
end