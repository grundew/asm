function [X, t] = convFreqWithTime(V, fs_v, x_ex, fs_ex)
%convFreqWithTime - convolves multiple frequency domain functions with a
%time domain function.
%
% [X, t] = convFreqWithTime(V, fs_v, x_ex, fs_ex)
% 
% Input:
% V     - One or more frequency domain signals
% fs_v  - Sampling frequency of coloumns V (scalar)
% x_ex  - Amplitude of time domain signal (vector)
% fs_ex - Sampling frequency of x_ex (scalar)
%
% Output:
% X - Convolved signal
% t - Time vector for which the coloumns of X are sampled


% Interpolate the excitation pulse to fs_v if sampling rate differs
if fs_v ~= fs_ex
    t_ex = (0:length(x_ex)-1)/fs_ex;
    t_v = 0:1/fs_v:max(t_ex);
    x_ex_ = interp1(t_ex, x_ex, t_v, 'linear');
else
    x_ex_ = x_ex;
end

nfft = size(V, 1);
X_ex_ = ifft(conj(hilbert(x_ex_)), nfft);
X_ex = repmat(X_ex_(:), 1, size(V, 2));
X = fft(V.*X_ex, nfft);
t = (0:nfft-1)/fs_v;
end