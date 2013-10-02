function P = planePistonPressureShiftedAngularSpectrum(kx, a, c, rho, deltax)
% P = planePistonPressureShiftedAngularSpectrum(kx, a, c, rho, deltax);
%   Calculates the angular spetrum of a plane piston shifted away from the z-axis.
%
% Input:
% kx - wave number (1/m)
% a - transducer radius (m)
% c - speed of sound (m/s)
% rho - density (kg/m3)
% deltax - shifting distance (m)

% Based on Ultrasonic Nondestructive Evaluation System - Lester W. Schmerr
% Jr and Sung-Jin Song.
v0 = 2/c/rho/2/pi/a^2;

% Transducer spacial sensitivities
P = 2*pi*besselj(0, kx*deltax).*planePistonPressureAngularSpectrum(kx, aTx, c, rho);
end
