function [x, t, T] = planeWaveTransmissionTimeSignal(model, xPulse, tPulse,...
    theta, dist, alphaLambda, nfft)
% [xT, t] = planeWaveTimeSignal(model, xPulse, tPulse, dist, alphaLambda)
%
% Calculates the transmission response of a solid plate embedded in a fluid
% and convolves with an excitation pulse.
%
% Input:
% model - MultiLayerModel object (Defines the material properties).
% xPulse - The amplitudes of the excitation pulse.
% tPulse - The time stamps of the excitation pulse samples.
% dist - Distance from target to receiver.
% alphaLambda - Absorption in the solid plate normalized to the wavelength.
% nfft - Number of points for the fft
%
% Output:
% xT - Transmitted time signal.
% t - Time.
% T - Transmission coefficient

if nargin < 4
    theta = 0;
    dist = 0;
    alphaLambda = 0;
    nfft = 2^14;
elseif nargin < 5
    dist = 0;
    alphaLambda = 0;
    nfft = 2^14;
elseif nargin < 6
    alphaLambda = 0;
    nfft = 2^14;
elseif nargin < 7
    nfft = 2^14;
end

fs = 1/(tPulse(2)-tPulse(1));
f = (0:nfft-1)*fs/nfft;

% Fourier transform of the pulse
Ypulse = ifft(xPulse, nfft);

% Get the reflection and transmission coefficients
alphaL = -alphaLambda*f/model.solid.v;
alphaL(f==0) = alphaL(2)*1e3;
T = transmissionCoefficientAnalytical(f, sin(theta), model, alphaL);

% Phase shift from propagating distance dist
kt = 2*pi*f/model.fluid(2).v;
Et = exp(1i*kt*dist);

% Convolve the pulse and the reflection/transmission coefficient
x = fft(Ypulse(:).*T(:).*Et(:), nfft);
t = (0:length(x)-1)/fs;
end