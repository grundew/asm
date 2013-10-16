function [xR, xT, t, R] = planeWaveTimeSignal(model, xPulse, tPulse, theta, dist, nfft)
% [xR, xT, t] = planeWaveTimeSignal(model, xPulse, tPulse)
%
% Calculates the response of a solid plate embedded in a fluid and
% convolves with an excitation pulse.
%
% Input:
% model - MultiLayerModel object (Defines the material properties).
% xPulse - The amplitudes of the excitation pulse.
% tPulse - The time stamps of the excitation pulse samples.
% theta - Angle of incidence.
% dist - Distance from target to receiver.
% nfft - Number of points for the fft
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
    nfft = 2^13;
elseif nargin < 5
    dist = 0;
    nfft = 2^13;
elseif nargin < 6
    nfft = 2^13;
end


%f = fftshift((0:nfft-1)*fs/nfft);
fs = 1/(tPulse(2)-tPulse(1));
df = fs/(nfft);
f = fftshift(((0:(nfft-1))'-floor(nfft/2))*df);

% Fourier transform of the pulse
Ypulse = ifft(xPulse, nfft);

% Get the reflection and transmission coefficients
R = fluidSolidFluidReflectionCoefficient(f, theta, model);

% Phase shift from propagating distance dist
kr = 2*pi*f/model.fluid(1).v;
kt = 2*pi*f/model.fluid(2).v;
Er = exp(1i*kr*dist);
Et = exp(1i*kt*dist);

% Convolve the pulse and the reflection/transmission coefficient
xR = fft(Ypulse(:).*R(:).*Er(:), nfft);
xT = fft(Ypulse(:).*R(:).*Et(:), nfft);
t = (0:length(xR)-1)/fs;
end