function [V, Wl, Wt, thetac_l, thetac_t] = fluidSolidReflectionCoefficient(theta, fluid, solid)
% Calculate reflection and transmission coeffiencts for a fluid-solid
% interface.
% 
% See Brekhovskikh, Acoustics of Layered Media I page 94.
% 
% Input:
% f: Frequncies (vector or scalar)
% theta: Incidence angle (vector or scalar)
% solid: Solid half-space properties
% fluid: Fluid half-space properties
% 
% The layer properties are structs with two fields, density and v
% (compatible with materials classes).
% 
% Output:
% V: Reflection coefficient
% Wl: Transmission coefficient, longitudenal
% Wt: Transmission coefficient, transverse
% thetac_l: Critical angle, longitudenal
% thetac_t: Critical angle, transversal
%
% Example:
% fluid = struct('v', 1500, 'density', 1000);
% solid = struct('v', 5349, 'vShear', 3158, 'density', 7810);
% f = 100e3:1e3:1000e3;
% theta = 0.01:0.001:0.2;
% [V, Wl, Wt, thetac_l, thetac_t] = fluidSolidReflectionCoefficient(f, theta, fluid, solid);

% Sound speeds
c = fluid.v;
cl1 = solid.v;
ct1 = solid.vShear;

% Densities
rho = fluid.density;
rho1 = solid.density;

% Angles in solid, l = longitudenal, t = transverse
theta_l = asin(cl1/c*sin(theta));
theta_t = asin(ct1/c*sin(theta));

% Impedances
Z = rho*c./cos(theta);
Zl = rho1*cl1./cos(theta_l);
Zt = rho1*ct1./cos(theta_t);

% Calculate reflection coefficient (Eq 4.2.24)
com = Zl.*cos(2*theta_t).^2 + Zt.*sin(2*theta_t).^2;
num = com - Z;
den = com + Z;
V = num./den;

% Transmission coefficient for longitudenal and shear wave in the solid.
% (Eq 4.2.25)
num = 2*rho/rho1*Zl.*cos(2*theta_t);
den = Zl.*cos(2*theta_t).^2 + Zt.*sin(2*theta_t).^2 + Z;
Wl = num./den;
% (Eq 4.2.25)
num = -2*rho/rho1.*Zt.*sin(2*theta_t);
Wt = num./den;

% Calculate critical angle for longitudenal and shear wave
thetac_l = asin(c/cl1);
thetac_t = asin(c/ct1);
end