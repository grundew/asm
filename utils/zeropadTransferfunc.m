function Y = zeropadTransferfunc(X)
% zeropadTransferFunc zero pads calculates the number of fft points and
% calculates the sampling frequency and frequency vector.
%
% Input:
% X - Matrix with transfer functions
% fx - Frequency vector with length = size(X, 1)
Y = [X; zeros(size(X, 1)-4, size(X, 2))];
end