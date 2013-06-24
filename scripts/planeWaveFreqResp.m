%% Plane wave normal incidence frequency response
fs = 6e6;
nfft = 2^15;
f = (0:nfft-1)*fs/nfft;

%% Excitation pulse
tend = 50e-6;
t = (0:1/fs:tend)';

% Start and stop frequencies
f0 = 400e3;
f1 = 1200e3;

% Window
wndw = rectwin(length(t));

% Real chirp
y = wndw.*chirp(t, f0, tend, f1, 'linear', 270);

% Pad with zeros
y = [zeros(100, 1); y; zeros(100, 1)];
t = (0:length(y)-1)/fs';

% Make the waveform Analytic
y = conj(hilbert(y));
Y = ifft(y, nfft);

%% Plot the pulse
figure
subplot(211)
plot(t, real(y))
subplot(212)
plot(f, db(abs(Y).^2))

%% Do the whole she bang
rho_fluid = 150;
v_fluid = 420;
v = 0;
d_Tx = 10e-3; % Doesn't matter
d_Rx = 10e-3; % Doesn't matter

for i = 1:length(v)
    %% Fluid-fluid-fluid model
    fluid1 = struct('v', v_fluid, 'density', rho_fluid);
    % fluid1 = struct('v', 340, 'density', 1000);
    fluid3 = fluid1;
    layer = struct('v', 5850, 'density', 7850, 'vShear', 3218);
    theta = v(i);
    
    d = 25e-3;
    fres = 0.5*5850/d;
    
    % Fluid-solid-fluid model
    model = MultiLayerModel(fluid1, layer, fluid3, d);
    
    [R, T] = analyticRTFast(f, theta, model);
    % T = transmissionCoefficientAnalytical(freq, q, model, alpha_L);
    
    %% Convolve the pulse and the reflection coefficient
    xt = fft(Y(:).*T(:), nfft);
    xr = fft(Y(:).*R(:), nfft);
    tx = (0:length(xt)-1)/fs;
    
    %% Plot both signals
    figure
    subplot(211)
    plot(tx, real(xt))
    legend('Transmitted')
    % plot(tx, real(xt), tx, real(xr));
    % legend('Transmitted', 'Reflected')
    title(sprintf('theta = %d', theta));
    ax = subplot(212);
    hold(ax, 'all');
    plotprd(xt, f/fres, ax);
    plotprd(xr, f/fres, ax);
    ylabel('dB')
    xlabel('Normalized frequency (f/f_1)')
end