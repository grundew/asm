function p = integrateTransmittedPressure(f, aRx, aTx, model, d1, d3)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
% f - Frequency (vector)
% aRx - Radius of plane piston receiver (scalar)
% aTx - Radius of plane piston transmitter (scalar)
% model - MultiLayerModel object
% R - Reflection or transmission coefficient
% d1 - Distance from transmitter to plate (scalar)
% d3 - Distance from receiver to plate (scalar)

w = 2*pi*f;
ntheta = 2^22;
thetamax = 0.5*pi;

% Front side fluid properties
c = model.fluid(1).v;
rho = model.fluid(1).density;

theta = linspace(0, thetamax, ntheta);
dtheta = theta(2) - theta(1);
q = sin(theta);

p = zeros(nfreq, 1);
for i = 1:nfreq
    kx = w(i)./c*q;
    kz = w(i)./c*sqrt(1-q.^2);
    
    % Transducer spacial sensitivities
    PhiRx = planePistonPressureAngularSpectrum(kx, aRx, c, rho);
    PhiTx = planePistonPressureAngularSpectrum(kx, aTx, c, rho);
    % Plate response, angular
    T = transmissionCoefficientAnalytical(f(i), q, model, alpha_L);
    % Phase shift from transmitter to plate and from plate to receiver
    Phase = exp(1i*kz*(d1 + d3));
    
    % Integrate over angles
    p(i) = integrate(T.*PhiRx.*PhiTx.*Phase, dtheta);
end

end

function p = integrate(I, h)
p = h*(0.5*(I(1) + I(end)) + sum(I(2:end-1)));
end