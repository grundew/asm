function [xR, xT, t, R, T] = planeWaveTimeSignal(model, xPulse, tPulse, theta, dist)
% [xR, xT, t] = planeWaveTimeSignal(model, xPulse, tPulse)
%
% Calculates the response of a solid plate embedded in a fluid and
% convolves with an excitation pulse.
%
% Input:
% model - MultiLayerModel object (Defines the material properties).
% xPulse - The amplitudes of the excitation pulse.
% tPulse - The time stamps of the excitation pulse samples.
% 
% Output:
% xR - Reflected time signal.
% xT - Transmitted time signal.
% t - Time.
% R - Reflection coefficient
% T - Transmission coefficient

if nargin < 4
    theta = 0;
    dist = 0;
elseif nargin < 5
    dist = 0;
end

fs = 1/(tPulse(2)-tPulse(1));
nfft = 2^13;
f = (0:nfft-1)*fs/nfft;

% Fourier transform of the pulse
Ypulse = ifft(xPulse, nfft);

% Get the reflection and transmission coefficients
[R, T] = fluidSolidFluidReflectionCoefficient(f, theta, model);

% Phase shift from propagating distance dist
kr = 2*pi*f/model.fluid(1).v;
kt = 2*pi*f/model.fluid(2).v;
Er = exp(1i*kr*dist);
Et = exp(1i*kt*dist);

% Convolve the pulse and the reflection/transmission coefficient
xR = fft(Ypulse(:).*R(:).*Er(:), nfft);
xT = fft(Ypulse(:).*T(:).*Et(:), nfft);
t = (0:length(xR)-1)/fs;
end