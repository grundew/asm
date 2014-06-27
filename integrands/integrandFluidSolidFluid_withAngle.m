function I = integrandFluidSolidFluid_withAngle(gamma, f, aRx, aTx,...
    c, rho, d, model, alpha)
% I = integrandFluidSolidFluid_withAngle(theta, f, aRx, aTx, c, rho, d, model, alpha)
%   orofinoIntegrand is the integrand in equation (28) in Ref. 1. This is
%   used as an input to the waveNumberIntegration function. It uses an
%   analytical expression for the transmission coefficient, assuming that
%   the fluid on both sides of the solid plate have the same properties.
% 
% Input:
% gamma - Angle (scalar or vector) between the plane wave and the
%         transducer surface normal.
% f - Frequency (scalar)
% aRx - Radius of the receiver
% aTx - Raidus of the transmitter
% c - Speed of sound in the propagation fluid
% rho - Density of the propagation fluid
% d - Distance from transmitter to the solid plate
% model - MultiLayerModel object
% alpha - The misalignment of the plate wrt a plane parallell to the
%         transducer plate.
% 
% Output:
% I - Integrand evaluated at theta
%
% References:
% 1. Orofino, 1992. http://dx.doi.org/10.1121/1.405408

% Angular frequency and total length of wave vector
c = model.fluid(1).v;
rho = model.fluid(1).density;

w = 2*pi*f;
k = w./c;

%% Transmitter spatial spectrum
q = sin(gamma);
p = sqrt(1-q.^2);
kx = k*q;
PhiTx = planePistonPressureAngularSpectrum(kx, aTx, c, rho);
idltgamma = gamma <= alpha;

%% Receiver spatial spectrum
% See confluence page on Angular Spectrum Method for details on the
% relations between the angles
gamma_rx = gamma - 2*alpha;
gamma_rx(idltgamma) = 2*alpha - gamma(idltgamma);
q_rx = sin(gamma_rx);
kx = k*q_rx;
PhiRx = planePistonPressureAngularSpectrum(kx, aRx, c, rho);

%% Plate response, angular
% See confluence page on Angular Spectrum Method for details on the
% relations between the angles
theta = gamma - alpha;
theta(idltgamma) = alpha - gamma(idltgamma);
R = analyticRTFast(f, theta, model);
% R = fluidSolidFluidReflectionCoefficient(f, theta, model);
% R = R.';

%% Phase shift from transmitter to plate and from plate to receiver
z = d + d*cos(2*alpha);
r = -d*sin(2*alpha);
kz = k*p;
kr = k*q;
Phase = exp(1i*(kr*r + kz*z));

I = R.*k.*q.*PhiRx.*PhiTx.*Phase.*k.*p;

if any(isnan(I))
    fprintf('NaN value detected at frequency %f and angle %f\n', f, gamma(isnan(I)));
end

end

