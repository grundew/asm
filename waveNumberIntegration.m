function [V, f] = waveNumberIntegration(model, integrandFunc, nfft, fs, thetamax)
% [V, f] = waveNumberIntegration(model, integrandFunc, nfft, fs, thetamax)
%
%   Function for integrating a spectrum of plane waves, see Equation
%   (29) in Ref 1.
%
% Input:
% model - MultiLayerModel object
% integrandFunc - Function handle to 
% nfft - Number of fft points
% fs - Sampling frequency
% thetamax - Limit of integration
%
% References:
% 1. Orofino, 1992. http://dx.doi.org/10.1121/1.405408

%% Samplings stuff
if nargin < 5
    thetamax = pi/2;
end
f = (0:nfft-1)*fs/nfft;

%% Transducer specs
aTx = model.a_Tx;
aRx = model.a_Rx;

%% Material parameters
d1 = model.d_Tx;
d3 = model.d_Rx;
rho = model.fluid(1).density;
c = model.fluid(1).v;
alphaLambda = model.alphaLambda;

%% Integrate over all angles for the point on the axis
nf = length(f);
V = zeros(nf, 1);

%% Loop over frequency and integrate over theta from zero to thetamax
for i = 1:nf
    func = @(theta) integrandFunc(theta, f(i),...
        aRx, aTx, c, rho, d1, d3, model, alphaLambda);
    V(i) = 2*pi*quadgk(func, 0, thetamax);
end

end