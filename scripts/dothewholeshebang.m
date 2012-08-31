%% Do the whole she bang!

%% Samplings stuff
fs = 2e6;
nfft = 2^15;
ntheta = 2^12;
thetamax = pi/2-0.01;
theta = linspace(-thetamax, thetamax, ntheta)';

% Observation point
z = 10e-2;
x = 0;

%% Transducer specs
a = 12e-3;

%% Material parameters
rho_fluid = 1000;
v_fluid = 1500;
v_layer = 5850;
fluid1 = struct('v', v_fluid, 'density', rho_fluid);
fluid3 = fluid1;
layer = struct('v', v_layer, 'density', 7850, 'vShear', 3218);
d = 10.15e-3;

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, d);
thresh = 1e-10;

%% Excitation pulse
tend = 50e-6;
f = fftshift((-nfft/2:nfft/2-1)*fs/nfft);
t = (0:1/fs:tend)';

% Start and stop frequencies
f0 = 200e3;
f1 = 500e3;

% Window
% wndw = rectwin(length(t));
wndw = gausswin(length(t));

% Real chirp
% y = wndw.*chirp(t, f0, tend, f1, 'linear', 270);

% Analytic chirp
alpha = (f1-f0)/tend;
y = wndw.*1i.*exp(-1j*2*pi*alpha*t.^2/2);

% Pad with zeros
y = [zeros(100, 1); y; zeros(100, 1)];
t = (0:length(y)-1)/fs';
Y = ifft(y, nfft);

% Plot the chirp with the spectrum
figure
subplot(211)
plot(t, real(y))
subplot(212)
plot(f, db(abs(Y)/max(abs(Y))), '.')

%% Integrate over all angles for the point on the axis
p = zeros(nfft, 1);
idfpos = find(f>0);
fpos = f(idfpos);
q = sin(theta);

% figure
for i = 1:length(fpos)
    freq = fpos(i);
    % Vertical and horizontal wave number
    kz = 2*pi*freq*cos(theta)/v_fluid;
    kx = 2*pi*freq*q/v_fluid;
    
    % Compute the factors in the integrand
    Phi = angularPlaneWaveSpectrumPiston(a, v_fluid, q, freq);
    [~, T] = fluidSolidFluidReflectionCoefficient(freq, theta, model, thresh);
    p(idfpos(i)) = propagateReflectedWave(x, z,...
        freq, q, kx, kz, rho_fluid, Phi, T);
    
    % Plot the integrands
    %     subplot(211)
    %     plot(theta, db(abs(T)/max(abs(T))))
    %     subplot(212)
    %     plot(theta, db(abs(Phi)/max(abs(Phi))), '.')
    %     pause(0.1)
end

%% Plot the frequency and impulse response of the observation point
figure
subplot(211)
plot(f, real(p), '.')
subplot(212)
plot(f, imag(p), '.')

figure
plot((0:nfft-1)/fs, real(fft(p, nfft)))

%% Convolve the pulse and the impulse response of the observation point
yy = conv(y, fft(p, nfft));
tt = (0:length(yy)-1)/fs;

%% Plot the convolved pulse
figure
plot(tt, real(yy))
xlabel('time')
ylabel('amplitude')
title('Transmitted pulse')