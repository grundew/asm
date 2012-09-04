%% Plane wave normal incidence frequency response
reflection = true;

%% Generate the chirp pulse
fs = 1e6;
tend = 50e-6;
t = (0:1/fs:tend)';
f0 = 100e3;
f1 = 500e3;

% Real chirp
% wndw = gausswin(length(t));
% y = wndw.*chirp(t, f0, tend, f1, 'linear', 270);

% Analytic chirp
alpha = (f1-f0)/tend;
wndw = gausswin(length(t));
y = wndw.*exp(-1j*2*pi*alpha*t.^2/2);
y = [zeros(100, 1); y; zeros(100, 1)]-mean(y);
t = (0:length(y)-1)/fs';

nfft = 2^20;
% f = fftshift((-nfft/2:nfft/2-1)*fs/nfft);
f = (0:nfft-1)*fs/nfft;
Y = ifft(y, nfft);

%% Plot the pulse
figure
subplot(211)
plot(t, real(y))
subplot(212)
plot(f, abs(Y))

%% Do the whole she bang
rho_fluid = 1000;
v_fluid = 1500;
v = 0.1;
for i = 1:length(v)
    %% Fluid-fluid-fluid model
    fluid1 = struct('v', v_fluid, 'density', rho_fluid);
    % fluid1 = struct('v', 340, 'density', 1000);
    fluid3 = fluid1;
    layer = struct('v', 5850, 'density', 7850, 'vShear', 3218);
    theta = v(i);
    %theta = 1e-4;
    d = 10.15e-3;
    
    % Fluid-solid-fluid model
    model = MultiLayerModel(fluid1, layer, fluid3, d);
    
    %     tic
    %     [V, W] = analyticRT(f, theta, model);
    %     toc
    tic
    [R, T] = fluidSolidFluidReflectionCoefficient(f, theta, model);
    toc
    
    %% Convolve the pulse and the reflection coefficient
    if reflection
        TT = R(:);
    else
        TT = T(:);
    end
    
    %     x = fft(Y(:).*TT(:), nfft);
    x = conv(y, fft(TT, nfft));
    xprd = abs(ifft(x.*hann(length(x)), nfft)).^2;
    tailstart = 800;
    tailend = 1200;
    xtailprd = abs(ifft(x(tailstart:tailend).*hann(tailend-tailstart+1), nfft)).^2;
    tx = (0:length(x)-1)/fs;
    
    if i == 1
        [axEcho, axTail] = plotSignal(tx, x, tx(tailstart:tailend), x(tailstart:tailend),...
            f, xprd, xtailprd, 0.5/d*[(1:2)*layer.v, 3*layer.vShear]);
    else    
        plotSignal(tx, x, tx(tailstart:end), x(tailstart:end),...
            f, xprd, xtailprd, 0.5/d*[(1:2)*layer.v, 3*layer.vShear], axEcho, axTail);
    end
     
end
