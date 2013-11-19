function I = integrandFluidSolidFluidTransmission_withLossAndDisplacement(k_r, w, aRx, aTx,...
    c, rho, d1, d3, model, x0, alphaLambda_dB)
% I = orofinoIntegrand(theta, f, aRx, aTx, c, rho, d1, d3, model, alphaLambda)
%   orofinoIntegrand is the integrand in equation (28) in Ref. 1. This is
%   used as an input to the waveNumberIntegration function. It uses an
%   analytical expression for the transmission coefficient, assuming that
%   the fluid on both sides of the solid plate have the same properties.
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
% model - MultiLayerModel object
% alphaLambda - Damping factor
% 
% Output:
% I - Integrand evaluated at theta
%
% References:
% 1. Orofino, 1992. http://dx.doi.org/10.1121/1.405408
%
k = w./c;

% k_r - input
k_z = sqrt(k^2 - k_r.^2);
sintheta_z = k_r/k;

%% Transducer spacial sensitivities
S_Rx = planePistonPressureAngularSpectrum(k_r, aRx, c, rho);
S_Tx = planePistonPressureAngularSpectrum(k_r, aTx, c, rho);

%% Displacement factor
dispRx = 2*pi*besselj(0, x0*kr);

%% Plate response, angular
% Multiply with wave length and convert from dB to linear
if alphaLambda_dB > 0
    alphaL = 10.^(alphaLambda_dB*w/2/pi/model.solid.v/20);
else
    alphaL = 0;
end

T = transmissionCoefficientAnalytical(w/2/pi, sintheta_z, model, alphaL);

%% Phase shift from transmitter to plate and from plate to receiver
Phase = exp(1i*k_z*(d1 + d3));
I = k_r.*dispRx.*k_z/k/rho/c.*Phase.*S_Rx.*S_Tx.*T;
end

