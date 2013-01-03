function [V, W] = reflectionCoefficientLayeredMedia(freq, theta, model)
% [V, W] = reflectionCoefficientLayeredMedia( input_args )
%

s = 1i*2*pi*freq;
q = sin(theta);

% Allocate parameters
nt = length(theta);
nf = length(freq);
V = zeros(nt, nf);
W = zeros(nt, nf);

% Speed of sounds
c_0 = model.fluid(1).v;
cShear_0 = model.fluid(1).vShear;

c_np1 = model.fluid(2).v;
cShear_np1 = modul.fluid(2).vShear;

vShear = model.solid.vShear;
vLong = model.solid.v;

% Densities
rho_0 = model.fluid(1).density;
rho_np1 = model.fluid(2).density;
rho_n = model.solid.density;

% Thickness of plate
d_n = model.thickness;

% Lamé coefficient
mu0 = rho_0*cShear_0^2;
mu_np1 = rho_np1*cShear_np1^2;
mu_n = rho_n*cShear_0^2;

for i = 1:nf
    % Calculate angle independent parameters here
    
    for j = 1:nt
        % Horizontal wave number component
        xi = s(i)*q(j)/c_0;
        
        % Calculate a_0
        k0 = sqrt(s(i)^2/c_0 - xi^2);
        gamma0 = rho_0*s^2 - 2*mu0*xi^2;
        a_0 = 2*rho_0*s^2*[gamma0/k0; -2*mu0*xi; xi/k0; -1];
        
        % Calculate a_{n+1}
        k_np1 = sqrt(s(i)^2/c_np1 - xi^2);
        h_np1 = sqrt(s(i)^2/cShear_np1 - xi^2);
        
        a_np1 = [k_np1, xi; -xi, h_np1; 2*mu_np1*xi, -gamma; -gamma, 2*mu_np1*xi*h_np1];
        
        % Calculate B_n
        % Make this a for-loop to take into account several layers
        v = 2*mu_n*xi^2/rho_n/s(i)^2;
        k_n = sqrt(s(i)^2/c_n - xi^2);
        h_n = sqrt(s(i)^2/cShear_n - xi^2);
        c_h = cosh(-h_n*d_n);
        c_k = cosh(-k_n*d_n);
        r_k = (k_n/xi)^2;
        r_h = (h_n/xi)^2;
        s_k = xi/k_n/sinh(-k_n*d);
        s_h = xi/h/sinh(-h*d);
        eta = rho_n*s(i)^2/xi;
        
        B11 = (1 - v)*c_k + v*c_n;
        B12 = v*s_k*r_k - (1 - v)*s_h;
        B13 = (c_k - c_h)/eta;
        B14 = (s_h + s_k*r_k)/eta;
        B21 = (1 - v)*s_k - v*s_h*r_h;
        B22 = (1 - v)*c_h + v*c_k;
        B23 = (s_k + s_h*r_h)/eta;
        
        g = a_0*B_layer*a_np1;

    end
end

