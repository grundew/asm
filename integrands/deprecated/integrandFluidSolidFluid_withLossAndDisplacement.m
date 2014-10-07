function I = integrandFluidSolidFluid_withLossAndDisplacement(theta_z, f, aRx, aTx,...
    c, d1, d3, model, x0, alphaLambda_dB, reflection)

%
% For some unknown reason this runs quicker than integrandeFluidSolidFluid_planepiston
%
% I = orofinoIntegrand(theta, f, aRx, aTx, c, rho, d1, d3, model, alphaLambda)
%   orofinoIntegrand is the integrand in equation (28) in Ref. 1. This is
%   used as an input to the waveNumberIntegration function. It uses an
%   analytical expression for the transmission coefficient, assuming that
%   the fluid on both sides of the solid plate have the same properties.
% 
% Input:
% theta_z - Wave number in the x,y-plane
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


%% Compute wave numbers
k = 2*pi*f./c;
k_z = k*cos(theta_z);
sintheta_z = sin(theta_z);
kr = k*sintheta_z;


%% Transmitter spatial sensitivity
% No angle adjustment
xx = kr*aTx;
W = besselj(1, xx)./xx;
W(xx==0) = 0.5;
PhiTx = 2*pi*aTx^2*W;


%% Receiver spatial sensitivity
xx = kr*aRx;
Wrx = besselj(1, xx)./xx;
Wrx(xx==0) = 0.5;
PhiRx = 2*pi*aRx^2*Wrx;


%% Displacement factor
if x0 > 0
    dispRx = besselj(0, x0*kr);
else
    dispRx = 1;
end

%% Plate response, angular
% Multiply with wave length and convert from dB to linear
if alphaLambda_dB > 0
    alphaL = 10.^(alphaLambda_dB*w/2/pi/model.solid.v/20);
else
    alphaL = 0;
end

if reflection
    Plate = analyticRTFast(f, theta_z, model);
else
    Plate = transmissionCoefficientAnalytical(f, sintheta_z, model, alphaL);
end


%% Phase shift from transmitter to plate and from plate to receiver
Phase = exp(1i*k_z*(d1 + d3));
I = kr.*dispRx.*k_z/k/c.*Phase.*PhiRx.*PhiTx.*Plate;


end
