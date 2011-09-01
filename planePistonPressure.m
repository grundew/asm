function [p, theta] = planePistonPressure(x, z, a, f, c, rho, U)
% Function calculating the pressure from a plane piston source at position
% x, z for a monochromatic wave with frequency f, propagating in a medium
% with propagation speed c and density \rho.
%
% The exp{i \omega t} is not included.
%
% See Fundamentals of acoustics, Kinsler et. al. (p. 182, eq, 7.4.17)
% Note that this is the far-field approximation.
%
% p = planePistonPressure(x, z, a, f, c, rho, U)
%
% x,z: Points where the field should be calculated. z is the aucoustical axis.
% a: Radius of piston source.
% f: Frequency of monochromatic wave.
% c: Propagation speed in medium.
% rho: Density of medium.
% U: Amplitude of harmonic motion of the plane piston.


if ~exist('U', 'var')
    U = 1;
end

k = 2*pi*f/c;
% Distance from source to field point.
r = sqrt(x.^2 + z.^2);
% Sine of angle from source to feild point.
sintheta = x./r;

% Axial part:
P_ax    = 1i*0.5*rho*c*U*a./r*k*a;
% Angular part:
vv = k*a*sintheta;
H_theta = 2*besselj(1, vv)./vv;

% Set the angularpart equal unity for theta equal zero
H_theta(sintheta==0) = 1;

p = P_ax.*H_theta.*exp(-1i*k*r);
theta = asind(sintheta);
end