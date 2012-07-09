function [R, T, KK, thetal, thetas, kL, kS, kzL, kzS] = analyticRT(freq, theta, model)

% L is half the thickness
L = 0.5*model.thickness;
rho_f = model.fluid.density;
rho_s = model.solid.density;
c_L = model.solid.v;
c_S = model.solid.vShear;
c_f = model.fluid.v;
  
T = zeros(length(freq), length(theta));
R = zeros(length(freq), length(theta));
KK = zeros(length(freq), length(theta));
kL = zeros(length(freq), 1);
kS = zeros(length(freq), 1);
kzL = zeros(length(freq), length(theta));
kzS = zeros(length(freq), length(theta));
thetal = zeros(length(freq), length(theta));
thetas = zeros(length(freq), length(theta));
for i = 1:length(freq);
    % Frequency
    f = freq(i);
    % Wave vector longitudenal
    h = 2*pi*f/c_L;
    % Wave vector shear
    k = 2*pi*f/c_S;
    % Debug variables
    kL(i) = h;
    kS(i) = k;
    for j = 1:length(theta)
        theta_f = theta(j);
        K = 2*pi*f/c_f*sin(theta_f);
        % Angle longitudenal
        theta_L = asin(K./h);
        % Angle shear
        theta_S = asin(K./k);
        
        % Debug variables
        KK(i, j) = K;
        thetal(i, j) = theta_L;
        thetas(i, j) = theta_S;
        kzL(i, j) = h*cos(theta_L);
        kzS(i, j) = k*cos(theta_S);
        
        a = rho_s*c_L*cos(theta_f)*cos(2*theta_S).^2;
        b = rho_s*c_S*cos(theta_f)*sin(2*theta_S).^2;
        N1 = a./(rho_f*c_f*cos(theta_L)*sin(2*h*L*cos(theta_L)));
        N2 = b./(rho_f*c_f*cos(theta_S)*sin(2*k*L*cos(theta_S)));
        N = N1 + N2;
        M1 = a./(rho_f*c_f*cos(theta_L)*tan(2*h*L*cos(theta_L)));
        M2 = b./(rho_f*c_f*cos(theta_S)*tan(2*k*L*cos(theta_S)));
        M = M1 + M2;

        T(i, j) = 2*N./(2*M + 1i*(M.^2 - N.^2 - 1));
        R(i, j) = 1i*(M.^2 - N.^2 + 1)./(2*M + 1i*(M.^2 - N.^2 - 1));
    end
    
end

T(freq==0, :) = ones(nnz(freq==0), length(theta));
R(freq==0, :) = zeros(nnz(freq==0), length(theta));
end