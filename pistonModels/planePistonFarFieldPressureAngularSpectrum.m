function [P, kx] = planePistonFarFieldPressureAngularSpectrum(x, z, f, a, c, rho, nfft)
% p = planePistonFarFieldPressureAngularSpectrum(x, z, f, a, c, rho, nfft)
%
% This function calculates the angular pressure spectrum from a plane
% piston in the far field.

% Omega
w = 2*pi*f;
% Wave number
k = w/c;
% Particle speed at transducer surface
v0 = 1;
theta = atan(x./z);
y = 0;
% Spatial sampling frequency
dx = x(2) - x(1);

% Pressure in x, z-domain
V = besselj(1, k*a*sin(theta))./(k*a*sin(theta));
idz = k*sin(theta)==0;
V(idz) = 0.5;
E = exp(1i*k*sqrt(x.^2 + y.^2 + z.^2))./sqrt(x.^2 + y.^2 + z.^2);
p = -1i*w*rho*a^2*v0*V.*E;

% Pressure in kx, kz-domain
if nargin < 7
    nfft = 2^15;
end

P = fft(p, nfft);
kx = (0:nfft-1)/dx/nfft;
end