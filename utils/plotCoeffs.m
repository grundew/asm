function plotCoeffs(x, param, d, X, titlestr, xlabelstr, varargin)
figure
ax(1) = subplot(211);
plot(x, real(X), varargin{:})
title(sprintf('%s - theta = %f, d = %f', titlestr, param, d))
xlabel(xlabelstr)
ylabel('Real')
ax(2) = subplot(212);
plot(x, imag(X), varargin{:})
xlabel(xlabelstr)
ylabel('Imag')
linkaxes(ax, 'x')
end