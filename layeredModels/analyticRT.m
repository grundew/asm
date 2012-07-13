function [R, T, N1, N2, M1, M2, alpha_L, alpha_S] = analyticRT(freq, theta, model)

% L is half the thickness
L = 0.5*model.thickness;
rho_f = model.fluid.density;
rho_s = model.solid.density;
c_L = model.solid.v;
c_S = model.solid.vShear;
c_f = model.fluid.v;

% Output parameters
nt = length(theta);
nf = length(freq);
T = zeros(nf, nt);
R = zeros(nf, nt);
N1 = zeros(nf, nt);
N2 = zeros(nf, nt);
M1 = zeros(nf, nt);
M2 = zeros(nf, nt);
alpha_L = zeros(nf, nt);
alpha_S = zeros(nf, nt);
for i = 1:length(freq);
    % Frequency
    f = freq(i);
    % Wave vector longitudenal
    h = 2*pi*f/c_L;
    % Wave vector shear
    k = 2*pi*f/c_S;

    for j = 1:length(theta)
        theta_f = theta(j);
        K = 2*pi*f/c_f*sin(theta_f);
        % Angle longitudenal
        theta_L = asin(K./h);
        % Angle shear
        theta_S = asin(K./k);
        
        a = rho_s*c_L*cos(theta_f)*cos(2*theta_S).^2;
        b = rho_s*c_S*cos(theta_f)*sin(2*theta_S).^2;
        N1(i, j) = a./(rho_f*c_f*cos(theta_L)*sin(2*h*L*cos(theta_L)));
        N2(i, j) = b./(rho_f*c_f*cos(theta_S)*sin(2*k*L*cos(theta_S)));
        N = N1(i, j) + N2(i, j);
        M1(i, j) = a./(rho_f*c_f*cos(theta_L)*tan(2*h*L*cos(theta_L)));
        M2(i, j) = b./(rho_f*c_f*cos(theta_S)*tan(2*k*L*cos(theta_S)));
        M = M1(i, j) + M2(i, j);
        
        % Debug parameters
        alpha_L(i, j) = h*L*cos(theta_L);
        alpha_S(i, j) = k*L*cos(theta_S);
        
        % Calculate reflection and transmission coefficients
        T(i, j) = 2*N./(2*M + 1i*(M.^2 - N.^2 - 1));
        R(i, j) = 1i*(M.^2 - N.^2 + 1)./(2*M + 1i*(M.^2 - N.^2 - 1));
    end
    
end

T(freq==0, :) = ones(nnz(freq==0), length(theta));
R(freq==0, :) = zeros(nnz(freq==0), length(theta));
end