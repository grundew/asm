function I = integrandPerfectReflector(theta_z, w, aRx, aTx,...
    c, rho, d1, d3, x0)
% I = integrandPerfectReflector(k_r, w, aRx, aTx,...
%    c, rho, d1, d3, x0)
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

q = sin(theta_z);
p = cos(theta_z);

kr = k*q;
kz = k*p;


%% Transmitter spatial sensitivity
% No angle adjustment
xx = kr*aTx;
W = besselj(1, xx)./xx;
W(xx==0) = 0.5;
PhiTx = 2*pi*aTx^2*W;

% Receiver spatial sensitivity
% No angle adjustment
xx = kr*aRx;
W = besselj(1, xx)./xx;
W(xx==0) = 0.5;
PhiRx = 2*pi*aRx^2*W;


%% Displacement factor
dispRx = 2*pi*besselj(0, x0*kr);


%% Phase shift from transmitter to plate and from plate to receiver
Phase = exp(1i*kz*(d1 + d3));
I = k.*q.*dispRx.*PhiRx.*PhiTx.*Phase.*k.*p.^2;
end