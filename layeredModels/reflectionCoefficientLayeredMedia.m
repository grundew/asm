function [Rp, Tp, Rs, Ts] = reflectionCoefficientLayeredMedia(freq, theta, model)
% [Rp, Tp, Rs, Ts] = reflectionCoefficientLayeredMedia(freq, theta, model)
%

s = 1i*2*pi*freq;
q = sin(theta);

% Allocate parameters
nt = length(theta);
nf = length(freq);
Rp = zeros(nt, nf);
Tp = zeros(nt, nf);
Rs = zeros(nt, nf);
Ts = zeros(nt, nf);

% Speed of sounds
c_0 = model.fluid(1).v;
cShear_0 = model.fluid(1).vShear;

c_np1 = model.fluid(2).v;
cShear_np1 = model.fluid(2).vShear;

cShear_n = model.solid.vShear;
c_n = model.solid.v;

% Densities
rho_0 = model.fluid(1).density;
rho_np1 = model.fluid(2).density;
rho_n = model.solid.density;

% Thickness of plate
d_n = model.thickness;

% Lamé coefficient
mu_0 = rho_0*cShear_0^2;
mu_np1 = rho_np1*cShear_np1^2;
mu_n = rho_n*cShear_0^2;

for i = 1:nf
    % Calculate angle independent parameters here
    
    for j = 1:nt
        % Horizontal wave number component
        xi = s(i)*q(j)/c_0;
        
        % Calculate a_0
        k_0 = sqrt(s(i)^2/c_0 - xi^2);
        h_0 = sqrt(s(i)^2/cShear_0 - xi^2); % Problem when cShear_0 is zero
        gamma_0 = rho_0*s(i)^2 - 2*mu_0*xi^2;
        A_0 = 2*rho_0*s(i)^2*[gamma_0/k_0, -2*mu_0*xi, xi/k_0, -1];
        
        % Calculate a_{n+1}
        k_np1 = sqrt(s(i)^2/c_np1 - xi^2);
        h_np1 = sqrt(s(i)^2/cShear_np1 - xi^2);
        gamma_np1 = rho_np1*s(i)^2 - 2*mu_np1*xi^2;
        A_np1 = [k_np1, xi; -xi, h_np1; 2*mu_np1*xi*k_np1, -gamma_np1; -gamma_np1, 2*mu_np1*xi*h_np1];
        
        %% Make this a for-loop to take into account several layers
        % Calculate all helper parameters
        v = 2*mu_n*xi^2/rho_n/s(i)^2;
        eta = rho_n*s(i)^2/xi;
        k_n = sqrt(s(i)^2/c_n - xi^2);
        h_n = sqrt(s(i)^2/cShear_n - xi^2);
        c_h = cosh(-h_n*d_n);
        c_k = cosh(-k_n*d_n);
        s_h = xi/h_n/sinh(-h_n*d_n);
        s_k = xi/k_n/sinh(-k_n*d_n);
        r_h = (h_n/xi)^2;
        r_k = (k_n/xi)^2;
        
        % Calculate <g>
        B(1, 1) = (1 - v)*c_k + v*c_h;
        B(1, 2) = v*s_k*r_k - (1 - v)*s_h;
        B(1, 3) = (c_k - c_h)/eta;
        B(1, 4) = (s_h + s_k*r_k)/eta;
        B(2, 1) = (1 - v)*s_k - v*s_h*r_h;
        B(2, 2) = (1 - v)*c_h + v*c_k;
        B(2, 3) = (s_k + s_h*r_h)/eta;
        B(2, 4) = eta*v*(1 - v)*(c_k - c_h);
        B(3, 1) = B(2, 4);
        B(3, 2) = eta*((1 - v)^2*s_h + v^2*s_k*r_k);
        B(3, 3) = B(2, 2);
        B(3, 4) = B(2, 1);
        B(4, 1) = eta*((1 - v)^2*s_k + v^2*s_h*r_h);
        B(4, 2) = B(1, 3);
        B(4, 3) = B(1, 2);
        B(4, 4) = B(1, 1);
        
        G = A_0*B*A_np1;
        
        % Calculate <j>
        C_0 = [4*mu_0^2*xi^2 - gamma_0^2/h_0/k_0, 4*mu_0*xi*gamma_0/k_0, 4*mu_0^2^xi^2 + gamma_0^2/h_0/k_0;...
            -rho_0*s(i)^2/k_0, 0, -rho_0*s(i)^2/k_0;...
            2*mu_0*xi + xi*gamma_0/h_0/k_0, 2*gamma_0/k_0, 2*mu_0*xi - xi*gamma_0/h_0/k_0;...
            2*mu_0*xi + xi*gamma_0/h_0/k_0, -4*mu_0*xi^2/k_0, 2*mu_0*xi - xi*gamma_0/h_0/k_0;...
            -rho_0*s(i)^2/h_0, 0, rho_0*s(i)^2/h_0;...
            -1 + xi^2/h_0/k_0, 2*xi/k_0, -(1 + xi^2/h_0/k_0)].';
        C_np1 = [h_np1*k_np1 + xi^2; -rho_np1*s(i)^2*k_np1;...
            xi*gamma_np1 - 2*mu_np1*xi*h_np1*k_np1; xi*gamma_np1 - 2*mu_np1*xi*h_np1*k_np1;...
            rho_np1*s(i)^2*h_np1; -4*mu_np1^2*xi^2*h_np1*k_np1 - gamma_np1^2];
        
        D(1, 1) = c_h*c_k + 2*v*(1 - v)*(1 - c_h*c_k) + ((1 - v)^2 + v^2*r_h*r_k)*s_h*s_k;
        D(1, 2) = (c_h*s_k + c_k*s_h*r_h)/eta;
        D(1, 3) = ((1 - 2*v)*(1 - c_h*c_k) - (1 - v - v*r_h*r_k)*s_h*s_k)/eta;
        D(1, 5) = -(c_h*s_k*r_k + c_k*s_h)/eta;
        D(1, 6) = (2*(1 - c_h*c_k) - s_h*s_k*(1 + r_h*r_k))/eta^2;
        D(2, 1) = eta*((1 - v)^2*c_k*s_h - v^2*c_h*s_k*r_k);
        D(2, 2) = c_h*c_k;
        D(2, 3) = v*c_h*s_k*r_k - (1 - v)*c_k*s_h;
        D(2, 5) = -s_h*s_k*r_k;
        D(3, 1) = eta*v*(1 - v)*(1 - 2*v)*(1 - c_h*c_k) + eta*((1 - v)^3 - v^3*r_h*r_k)*s_h*s_k;
        D(3, 2) = (1 - v)*c_h*s_k - v*c_k*s_h*r_h;
        D(5, 1) = -eta*(v^2*c_k*s_h*r_h + (1 - v)^2*c_h*s_k);
        D(6, 1) = eta^2*(2*v^2*(1 - v)^2*(1 - c_h*c_k) - (v^4*r_h*r_k + (1 - v)^4)*s_h*s_k);
        D(5, 2) = -s_h*s_k*r_h;
        D(1, 4) = D(1, 3);
        D(2, 4) = D(2, 3);
        D(2, 6) = D(5, 1);
        D(3, 4) = D(2, 2) - D(1, 1);
        D(3, 3) = D(3, 4) + 1;
        D(3, 5) = D(4, 2);
        D(4, 1) = D(3, 1);
        D(3, 6) = D(4, 1);
        D(4, 2) = D(3, 2);
        D(4, 3) = D(3, 4);
        D(4, 4) = D(3, 3);
        D(4, 5) = D(3, 2);
        D(4, 6) = D(3, 1);
        D(5, 3) = D(2, 4);
        D(5, 4) = D(2, 3);
        D(5, 5) = D(2, 2);
        D(5, 6) = D(2, 1);
        D(6, 2) = D(1, 5);
        D(6, 3) = D(1, 4);
        D(6, 4) = D(1, 3);
        D(6, 5) = D(1, 2);
        D(6, 6) = D(1, 1);
        
        J = C_0*D*C_np1;
        
        % Calculate reflection and transmission coefficient
        Rp(i, j) = J(1)/J(3);
        Rs(i, j) = J(2)/J(3);
        Tp(i, j) = G(1)/J(3);
        Ts(i, j) = G(2)/J(3);
    end
end

