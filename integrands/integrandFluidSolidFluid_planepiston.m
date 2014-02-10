function I = integrandFluidSolidFluid_planepiston(theta, f, aRx, aTx,...
    c, d1, d3, model, x0, alphaLambda_dB, reflection)
% I = orofinoIntegrand(theta, f, aRx, aTx,...
%                      c, rho, d1, d3, model,...
%                      alphaLambda_dB, reflection, alpha_plate)
%   orofinoIntegrand is the integrand in equation (28) in Ref. 1. This is
%   used as an input to the waveNumberIntegration function. It uses an
%   analytical expression for the transmission coefficient, assuming that
%   the fluid on both sides of the solid plate have the same properties.
% 
% Input:
% theta - Angle (scalar or vector)
% f - Frequency (scalar)
% aRx - Radius of the receiver
% aTx - Raidus of the transmitter
% c - Speed of sound in the propagation fluid
% rho - Density of the propagation fluid
% d1 - Distance from transmitter to the solid plate
% d3 - Distance from receiver to the solid plate
% x0 - Displacement of receiver from the acoustical axis of the transmitter
% model - MultiLayerModel object
% alphaLambda - Damping factor
% reflection - (Boolean) True if reflection coeffecient should be used
% alpha_plate - Angle of misalignment of the plate
%
% Output:
% I - Integrand evaluated at theta
%
% References:
% 1. Orofino, 1992. http://dx.doi.org/10.1121/1.405408

q = sin(theta);
p = sqrt(1-q.^2);
w = 2*pi*f;
k = w./c;
kr = k*q;
kz = k*p;

%% Transmitter spatial sensitivity
% No angle adjustment
xx = kr*aTx;
W = besselj(1, xx)./xx;
W(xx==0) = 0.5;
PhiTx = 2*W;

%% Receiver spatial sensitivities
xx = kr*aRx;
W = besselj(1, xx)./xx;
W(xx==0) = 0.5;
PhiRx = 2*W;

%% Plate response, angular
% Multiply with wave length and convert from dB to linear
%% Loss parameter
if alphaLambda_dB > 0
    alphaL = 10.^(alphaLambda_dB*f/model.solid.v/20);
else
    alphaL = 0;
end

%% Displacement factor
if x0 > 0
    dispRx = 2*pi*besselj(0, x0*kr);
else
    dispRx = 2*pi;
end

%% Reflection/Transmission coefficient
if reflection
    % TODO: Add loss
    Plate = analyticRTFast(w/2/pi, theta, model);
else
    Plate = transmissionCoefficientAnalytical(f, q, model, alphaL);
end

%% Phase shift from transmitter to plate and from plate to receiver
Phase = exp(1i*kz*(d1 + d3));

%% Assemble integrand
I = Plate.*k.*q.*dispRx.*PhiRx.*PhiTx.*Phase.*k.*p;
end