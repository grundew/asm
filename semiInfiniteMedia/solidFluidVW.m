function [Vll, Vlt, Vtt, Vtl, Wll, Wtl, theta, theta_t] = solidFluidVW(theta_l, rho_fluid, c_fluid, rho_solid, c_l, c_t)
% [Vll, Vlt, Vtt, Vtl, Wll, Wtl] = solidFluidVW(theta_l, rho_fluid, c_fluid, rho_solid, c_l, c_t)
% 
% Calculates the reflection and transmission of a l and t-wave of a sound
% wave in a fluid on a solid. See chapter 4 in Brekhovskikh.

% Angles
theta = asin(c_fluid/c_l.*sin(theta_l));
theta_t = asin(c_t/c_l.*sin(theta_l));

% Impedances
Z = rho_fluid*c_fluid./cos(theta);
Zl = rho_solid*c_l./cos(theta_l);
Zt = rho_solid*c_t./cos(theta_t);

% Helper variables
cos2t_t = cos(2.*theta_t);
sin2t_t = sin(2.*theta_t);
sint_t = sin(theta_t);
denom = (Z + Zt.*sin2t_t.^2 + Zl.*cos2t_t.^2);

% Reflection and transmission coefficients
Vll = (Z + Zt.*sin2t_t.^2 - Zl.*cos2t_t.^2)./denom;
Vlt = -(2.*(1-Vll).*cot(theta_l).*sint_t.^2)./cos2t_t;
Vtt = -(Z + Zl.*cos2t_t.^2 - Zt.*sin2t_t.^2)./denom;
Vtl = (1 + Vtt).*tan(theta_l).*cos2t_t/2./sint_t.^2;
Wll = (1 - Vll).*tan(theta).*cot(theta_l)./cos2t_t;
Wtl = (1 + Vtt).*tan(theta)/2./sint_t.^2;
end