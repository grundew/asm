function [X, kx, ky] = circularAperture(r, verbose)

%% Circular aperture
%r = 0.06; % Meter
dx = 0.001; % Meter
n = 1000;
x = zeros(n, n);
for i = 0:n-1
    for j = 0:n-1
        if sqrt(i*i*dx*dx + j*j*dx*dx) < r
            x(i+1, j+1) = 1;
        end
    end
end

x = [rot90(x,2) flipud(x); fliplr(x) x];

%% Calculate the FFT and the k-vector
nfft = 2^12;
X = fft2(x, nfft, nfft);
Xf = fftshift(X);
p = -nfft/2:nfft/2-1;
ky = p*(2*pi/dx/nfft);
kx = p*(2*pi/dx/nfft);

if exist('verbose', 'var') && verbose == true
    % Plot aperture
    figure;
    imagesc(x)
    axis equal
    
    
    %% Plot FFT in linear and logarithmic colorscale
    figure;
    ax(1) = subplot(211);
    axis xy
    axis equal
    imagesc(kx, ky, abs(Xf));
    axis xy
    axis equal
    ax(2) = subplot(212);
    imagesc(kx, ky, 20*log10(abs(Xf)));
    axis xy
    axis equal
    linkaxes(ax, 'xy')
    
    % Highligth the first side lobe
    theta = linspace(0, 2*pi, n);
    xsl = 5.14/r*cos(theta);
    ysl = 5.14/r*sin(theta);
    hold(ax(1), 'all')
    plot(ax(1), xsl, ysl, 'r.')
    hold(ax(2), 'all')
    plot(ax(2), xsl, ysl, 'b.')
    
    %% Plot the centre line
    figure
    id = size(X, 2)/2;
    plot(kx, 20*log10(abs(Xf(id, :))/max(abs(Xf(id, :)))), '.')
end