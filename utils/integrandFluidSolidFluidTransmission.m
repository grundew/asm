function I = integrandFluidSolidFluidTransmission(theta, f, aRx, aTx,...
    c, rho, d1, d3, model)
% I = integrandFluidSolidFluidTransmission(theta, f, aRx, aTx, c, rho, d1, d3, model, alphaLambda)
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
% model - MultiLayerModel object
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
kx = k*q;
kz = k*p;

%% Transducer spacial sensitivities
PhiRx = planePistonPressureAngularSpectrum(kx, aRx, c, rho);
PhiTx = planePistonPressureAngularSpectrum(kx, aTx, c, rho);

%% Plate response, angular
[~, T] = analyticRTFast(freq, theta, model);

%% Phase shift from transmitter to plate and from plate to receiver
Phase = exp(1i*kz*(d1 + d3));

I = T.*k.*q.*PhiRx.*PhiTx.*Phase.*k.*p;
end

