function ax = plotIR(f, R, ir, ax)
% Function for plotting real and imaginary part of a reflection or
% transmission coefficient.
% 
% plotR(f, R, ir, ax)
%
% f - Frequencies
% R - Reflection or transmission coefficient
% ir - If true, the real and imaginary part is plotted (default), otherwise
% the absolute value and the angle is plotted.
% ax - Vector of length two of axis handles (Optional). Axis is held before plotting.

% Check if we have some axes in the input.
if nargin == 4
    ish(1) = ishold(ax(1));
    ish(2) = ishold(ax(2));
    hold(ax(1), 'all')
    hold(ax(2), 'all')
else
    % If not create some axes.
    fig = figure();
    ax(1) = subplot(2, 1, 1, 'Parent', fig);
    ax(2) = subplot(2, 1, 2, 'Parent', fig);
    ish(1) = 0;
    ish(2) = 0;
end

% Check if user wants the real/imaginary or abs/phase.
if ~exist('ir', 'var') || isempty(ir) || ir
    plot(ax(1), f, real(R))
    plot(ax(2), f, imag(R))
    ylabel(ax(1), 'Real')
    ylabel(ax(2), 'Imag')
else
    plot(ax(1), f, abs(R))
    plot(ax(2), f, unwrap(angle(R)))
    ylabel(ax(1), 'Abs')
    ylabel(ax(2), 'Angle')
end

xlabel(ax(2), 'Frequency')

% Hold off the axes again if they weren't.
if ~ish(1)
    hold(ax(1), 'off')
end

if ~ish(2)
    hold(ax(2), 'off')
end

linkaxes(ax, 'x')
end
