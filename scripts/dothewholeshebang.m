%% Do the whole she bang!
debug = false;
saveresults = true;

%% Samplings stuff
fs = 2e6;
nfft = 2^12;
nq = 2^15;
qmax = 0.3;
q = linspace(-qmax, qmax, nq)';
theta = asin(q);
f = (0:nfft-1)*fs/nfft;
% f = 2.05*fres;
% f = [100e3, 1000e3, 2000e3];

%% Transducer specs
a = 9e-3;
arx = 3e-3;

% Observation point
xmax = 10e-2;
nx = 2^6;
h = 10e-2;
z = h;
xo = linspace(0, arx, nx);
% xo = 0;
zo = h;

%% Material parameters
rho_fluid = 1.5;
v_fluid = 350;
v_layer = 5850;
fluid1 = struct('v', v_fluid, 'density', rho_fluid);
fluid3 = fluid1;
layer = struct('v', v_layer, 'density', 7850, 'vShear', 3218);
d = 10e-3;
fres = 0.5*v_layer/10.15e-3;

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, d);
thresh = 1e-10;

%% Excitation pulse
tend = 50e-6;
t = (0:1/fs:tend)';

% Start and stop frequencies
f0 = 200e3;
f1 = 800e3;

% Window
% wndw = rectwin(length(t));
wndw = gausswin(length(t));

% Real chirp
y = wndw.*chirp(t, f0, tend, f1, 'linear', 270);

% Analytic chirp
% alpha = (f1-f0)/tend;
% y = wndw.*exp(-1j*2*pi*alpha*t.^2/2);
y = conj(hilbert(y));

% Pad with zeros
y = [zeros(100, 1); y; zeros(100, 1)];
t = (0:length(y)-1)/fs';
Y = ifft(y, nfft);

%% Integrate over all angles for the point on the axis
% q = sin(theta);

tic
% figure
nf = length(f);
pt = zeros(nf, nx);
pr = zeros(nf, nx);
for i = 1:nf
    % Time it
    if i == 1
        datestr(now)
        tic
    end
    
    freq = f(i);
    k = 2*pi*freq/v_fluid;
    
    % Vertical and horizontal wave number
    kz = k*cos(theta)/v_fluid;
    kx = k*q/v_fluid;
    
    % Compute the factors in the integrand
    Phi = planePistonPressureAngularSpectrum(z, kx, kz, a, v_fluid, rho_fluid);
    
    % Reflection or transmission coefficient
    % [R, T] = fluidSolidFluidReflectionCoefficient(freq, theta, model, thresh);
    [R, T] = analyticRTFast(freq, theta, model);

    Ht = Phi.*T.*exp(1i*kz*zo);
    Hr = Phi.*R.*exp(1i*kz*zo);

    for j = 1:nx
        E = exp(1i*kx.*xo(j));
        pt(i, j) = trapz(q, Ht.*E, 1);
        pr(i, j) = trapz(q, Hr.*E, 1);
    end
    
    % Time it
    if i == 2
        tme = 0.5*toc/60*length(f);
        fprintf('Estimated time of arrival: %f min\n', tme)
    end
    
    if debug
        debugplots(theta, freq, Phi, T, Ht, E)
    end
end

% % Integrate over position
pint = zeros(nfft, 1);
for i = 1:nfft
    pint(i) = trapz(xo, 2*pi*xo.*pt(i, :));
end
 
%% Convolve the pulse and the impulse response of the observation point
yt = conv(real(y), fft(pint, nfft));
yr = conv(real(y), fft(pint, nfft));
tt = (0:length(yt)-1)/fs;

%% Finnished
if saveresults
    dtstr = datestr(now, 'dd-mm-yyyy_HH-MM-SS');
    outfile = sprintf('wholeshebang_%s_thickness%d.mat', dtstr, d);
    fprintf('Finnished in %s min - %s\n', toc/60, dtstr)
    fprintf('Saved to %s\n', outfile)
    save(outfile);
end