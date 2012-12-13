function [V, Wl, Wt, theta_l, theta_t] = fluidSolidVW(theta, rho_fluid, c_fluid, rho_solid, c_l, c_t)
% [V, Wl, Wl] = fluidSolidRT(f, theta, rho_fluid, c_fluid, rho_solid, cl, ct)
% 
% Calculates the reflection and transmission of a l and t-wave of a sound
% wave in a fluid on a solid. See chapter 4 in Brekhovskikh.

% Angles
theta_l = asin(c_l/c_fluid.*sin(theta));
theta_t = asin(c_t/c_fluid.*sin(theta));

% Impedances
Z = rho_fluid*c_fluid./cos(theta);
Zl = rho_solid*c_l./cos(theta_l);
Zt = rho_solid*c_t./cos(theta_t);

% Helper variables
cos2t_t = cos(2.*theta_t);
sin2t_t = sin(2.*theta_t);
denom = (Zl.*cos2t_t.^2 + Zt.*sin2t_t.^2 + Z);

% Reflection and transmission coefficients
V = (Zl.*cos2t_t.^2 + Zt.*sin2t_t.^2 - Z)./denom;
Wl = (2*rho_fluid/rho_solid*Zl.*cos2t_t)./denom;
Wt = (-2*rho_fluid/rho_solid*Zt.*sin2t_t)./denom;
end