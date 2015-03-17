function I = integrandPerfectReflector(theta_z, w, aRx, aTx,...
    c, rho, d1, d3, ~, x0, ~)
% I = integrandPerfectReflector(k_r, w, aRx, aTx,...
%    c, rho, d1, d3, ~, x0, ~)
%
% This function returns the integrand for a pitch-catch setup on a perfect
% reflector, with possibly displaced receiver.
% 
% Input:
% k_r - Wave number in the x,y-plane
% w - Angular frequency (scalar)
% aRx - Radius of the receiver
% aTx - Raidus of the transmitter
% c - Speed of sound in the propagation fluid
% rho - Density of the propagation fluid
% d1 - Distance from transmitter to the solid plate
% d3 - Distance from receiver to the solid plate
% model - Not in use, but kept in input arguments to keep it compatible
%         with startAsmSimulation.
% alphaLambda - Not in use, but kept in input arguments to keep it compatible
%         with startAsmSimulation.
% 
% Output:
% I - Integrand evaluated at theta
%
k = w./c;

k_r = k*sin(theta_z);
k_z = k*cos(theta_z);

%% Transducer spacial sensitivities
S_Rx = planePistonPressureAngularSpectrum(k_r, aRx, c, rho);
S_Tx = planePistonPressureAngularSpectrum(k_r, aTx, c, rho);

%% Displacement factor
dispRx = 2*pi*besselj(0, x0*k_r);

%% Phase shift from transmitter to plate and from plate to receiver
Phase = exp(1i*k_z*(d1 + d3));
I = k_r.*dispRx.*k_z/k/rho/c.*Phase.*S_Rx.*S_Tx;
end