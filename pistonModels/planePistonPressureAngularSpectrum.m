function P = planePistonPressureAngularSpectrum(z, kx, kz, a, c, rho)
% p = planePistonPressureAngularSpectrum(theta, z, f, a, c, rho)
%
% Based on Ultrasonic Nondestructive Evaluation System - Lester W. Schmerr
% Jr and Sung-Jin Song.

% Normalize to one for numerical concerns
v0 = 2/c/rho/2/pi/a^2;

% Equation 8.10 and equation just after 8.6 from ref above
x = kx*a;
W = besselj(1, x)./x;
W(x==0) = 0.5;
V = c*rho*2*pi*a^2*v0*W;

% Phase factor
phase = exp(1i*kz.*z);
P = V.*phase;
end