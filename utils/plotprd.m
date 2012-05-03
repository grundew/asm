function ax = plotprd(x, tailstart, len, nfft, fs, ax)

x = x(:);
xtail = x(tailstart:tailstart+len-1);
% wndw = hann(length(xtail));
wndw = rectwin(length(xtail));
Xtail = abs(fft(xtail.*wndw, nfft)).^2;
f = (-nfft/2:nfft/2-1)*(fs/nfft);

if nargin < 6
    figure;
    ax = axes;
end
plot(ax, f, Xtail.^2);

end