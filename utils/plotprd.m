function ax = plotprd(x, f, ax)

if nargin < 3
    figure
    ax = axes();
end

nfft = length(f);
X = 1/nfft*abs(ifft(x, nfft)).^2;
plot(ax, f, db(X));

end