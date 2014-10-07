function [V, W, debug] =  fluidLayerReflectionCoefficient(f, theta, fluid3, fluid2, fluid1, d)
% Calculates the reflection and transmission coefficient for a fluid layer
% embedded in a fluid (same fluid on each side of the fluid layer).
%
% See Brekhovskikh, Acoustics of Layered Media I page 28.
% 
% Medium 3 is the medium with the incoming and reflected wave.
% Medium 2 is the layer with thickness d.
% Medium 1 is the medium with the transmitted wave.
%
% Input:
% f: Frequncies (vector or scalar)
% theta: Incidence angle (vector or scalar)
% fluid3: Medium 3 properties
% fluid2: Layer properties
% fluid1: Medium 1 properties
% 
% The layer properties are structs with two fields, density and v
% (compatible with materials classes).
%
% Example:
% fluid1 = struct('v', 1500, 'density', 1000);
% fluid3 = fluid1;
% layer = struct('v', 5900, 'density', 7850);
% f = 100e3:1e3:1000e3;
% theta = 0.01:0.001:0.2;
% d = 12.4e-3;
% [V, W] = fluidLayerReflectionCoefficient(f, theta, fluid3, layer, fluid1, d);
% plot(f, real(V(:, 1))
% hold all
% plot(f, real(W(:, 1))

f = f(:);
theta = theta(:).';
if numel(theta) > 1
    w = 2*pi*repmat(f, 1, length(theta));
else
    w = 2*pi*f;
end

if numel(f) > 1
    theta = repmat(theta, length(f), 1);
end
rho1 = fluid1.density;
rho2 = fluid2.density;
rho3 = fluid3.density;

c1 = fluid1.v;
c2 = fluid2.v;
c3 = fluid3.v;

% Compute length of wavevector in the mediums and horizontal wavenumber to
% get the angles
k2 = w/c2;
k3 = w/c3;
K = k3.*sin(theta);

theta2 = asin(K./k2);
theta3 = asin(K./k3);

% Calculate the impedances of the layers
Z1 = rho1*c1./cos(theta); 
Z2 = rho2*c2./cos(theta2);
Z3 = rho3*c3./cos(theta3);
fi = k2*d.*cos(theta2);

% Reflection and transmission coefficients (equation 2.4.9 and 2.4.13 page
% 28 in Brekhovskikh)
num = (Z1 + Z2).*(Z2 - Z3).*exp(-1i*2*fi) + (Z1 - Z2).*(Z2 + Z3);
den = (Z1 + Z2).*(Z2 + Z3).*exp(-1i*2*fi) + (Z1 - Z2).*(Z2 - Z3);
V = num./den;

W = (1 + V)./(cos(fi) - 1i*Z2.*sin(fi)./Z1);


% If f = 0 is present, set V = 0 and T = 1 at f = 0
V(f==0) = 0;
W(f==0) = 1;
end