function [RR] = angleToWaveNumber(f, theta, R, k_hor)
% Rearranges a reflection coefficient matrix sampled in f and theta to f
% and k_horizontal.
%
% Input:
% f: (vector) Frequencies for which the R is calculated.
% theta: (vector) Angles for which the R is calculated.
% K: (matrix) Horizontal wavenumber for which the R is calculated.
% R: (matrix) Reflection or transmission coefficient matrix.
%
% Output:
% RR: (matrix) Rearranged reflection or transmission matrix.
% K_hor: (vector) Horizontal wavenumbers for which the R is now sampled.

nf = length(f);
nk = length(k_hor);

RR = zeros(nf, nk);

% fig = figure;
% ax1 = subplot(211);
% ax2 = subplot(212);

for i = 1:nf
    k = 2*pi*f(i)/1500*sin(theta);
    RR(i, :) = interp1(k, R(i, :), k_hor);
    
    % Plot the original
    % plot(ax1, k, real(R(i, :)));
    % plot(ax2, k, imag(R(i, :)))
    % hold(ax1, 'all')
    % hold(ax2, 'all')

    % Plot the interpolated
    % plot(ax1, k_hor, real(RR(i, :)));
    % plot(ax2, k_hor, imag(RR(i, :)));
    % hold(ax1, 'off')
    % hold(ax2, 'off')
    
    % title(ax1, sprintf('Frequency: %d', f(i)));
    % Save
    % filename = fullfile('figures', sprintf('r%d.png', f(i)));
    %saveas(fig, filename, 'png')
end