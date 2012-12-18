function I = orofiniIntegrand(theta, f, aRx, aTx,...
    c, rho, d1, d3, model, alphaLambda)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
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
alphaL = -alphaLambda*f/model.solid.v;
T = transmissionCoefficientAnalytical(f, q, model, alphaL);

%% Phase shift from transmitter to plate and from plate to receiver
Phase = exp(1i*kz*(d1 + d3));

I = T.*k.*q.*PhiRx.*PhiTx.*Phase.*k.*p;
end

