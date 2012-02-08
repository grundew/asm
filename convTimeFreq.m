function [c, t] = convTimeFreq(x, fs, nfft, Y, fy, complex)
% c = convTimeFreq(x, fs, nfft, Y, f)
%
% Calculates the convolution between a time-domain signal (x) and a frequency-domain signal (Y)
% using fft on x, and then ifft on X*Y.
%
% c - Convoluted signal
% x - Time-domain signal
% fs - Sampling frequency for x
% nfft - Number of fft samples
% Y - Frequency-domain signal
% fy - Frequencies where Y is computed
% complex - Set to true if the imaginary part should be kept.

X = fft(x, nfft);
fx = fftshift((-nfft/2:nfft/2-1)*fs/nfft);
% Interpolate Y's to where X is sampled
if (length(fx) ~= length(fy)) || (length(fx) == length(fy) && sum(nnz(fx-fy)) > 0)
    Yx = interp1(fy, Y, fx, 'cubic', 0).';
else
    Yx = Y;
end

nifft = length(X);
c = ifft(X.*Yx, nifft);

t = (0:nfft-1)/fs;

if exist('complex', 'var') && ~complex
    c = real(c);
end

end