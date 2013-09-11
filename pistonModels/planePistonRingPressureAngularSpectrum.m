function P = planePistonRingPressureAngularSpectrum(kx, r_min, r_max, c, rho)
% p = planePistonPressureAngularSpectrum(theta, z, f, r_min, r_max, c, rho)
%
% Based on Ultrasonic Nondestructive Evaluation System - Lester W. Schmerr
% Jr and Sung-Jin Song.

% Normalize to one for numerical concerns
v0 = 2/c/rho/2/pi;

% Equation 8.10 and equation just after 8.6 from ref above
x1 = kx*r_max;
W1 = r_max*besselj(1, x1);
x2 = kx*r_min;
W2 = r_min*besselj(1, x2);
W1(x1==0) = 0.5;
W2(x2==0) = 0.5;
P = c*rho*2*pi*v0*(W1 - W2)./kx;
end