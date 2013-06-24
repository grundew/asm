%% Do the whole she bang!
debug = false;
saveresults = true;

%% Samplings stuff
fs = 2e6;
nfft = 2^4;
thetamax = pi/2;
f = (0:nfft-1)*fs/nfft;

%% Transducer specs
aTx = 5e-3;
aRx = 5e-3;
d1 = 10e-2;
d3 = 10e-2;

%% Material parameters
% Water layer
rho_fluid = 1000;
v_fluid = 1500;
fluid1 = struct('cp', v_fluid, 'rho', rho_fluid);

% Steel layer
v_layer = 5850;
vs_layer = 3218;
rho_layer = 7850;
poisson1 = (2*vs_layer^2 - v_layer^2)/2/(vs_layer^2 - v_layer^2);
d = 12e-3;
fres = 0.5*v_layer/d;
layer1 = struct('cp', v_layer, 'rho', rho_layer, 'poisson', poisson1, 'thickness', d);

% Concrete layer
cp2 = 3230;
cs2 = 1500;
poisson2 = (2*cs2^2 - cp2^2)/2/(cs2^2 - cp2^2);
layer2 = struct('cp', cp2, 'rho', 2600, 'poisson', poisson2,'thickness', 12e-3); % between fast and slow cement, thickness of last layer is ignored

% Model class

% Water - steel - concrete
% fluid = fluid1;
% layer = [layer1, layer2];
% msmpc = MultiShearModelPlateClass2(f(1), 1e-4, layer, fluid); % Initiliaze some random angle

% Water - steel - water
fluid = [fluid1, fluid1];
layer = layer1;
msmpc = MultiShearModelPlateClass(0, 0, layer, fluid); % Initiliaze some random angle
msmpc.computeShearVelocity();

%% Excitation pulse
tend = 50e-6;
t = (0:1/fs:tend)';

% Start and stop frequencies
f0 = 400e3;
f1 = 1200e3;

% Window
wndw = rectwin(length(t));
% wndw = gausswin(length(t));

% Real chirp
y = wndw.*chirp(t, f0, tend, f1, 'linear', 270);

% Analytic chirp
% alpha = (f1-f0)/tend;
% y = wndw.*exp(-1j*2*pi*alpha*t.^2/2);
y = conj(hilbert(y));

% Pad with zeros
t = (0:length(y)-1)/fs';
Y = ifft(y, nfft);

%% Integrate over all angles for the point on the axis
tic
% figure
nf = length(f);
pt = zeros(nf, 1);

for i = 1:nf
    % Time it
    if i == 1
        fprintf('Started: %s\n', datestr(now, 'dd-mm-yyyy_HH-MM-SS'));
        tic
    end
    
    freq = f(i);
    
    fun = @(xx) multiShearModelPlateClassIntegrand(xx, freq, aRx, aTx,...
        v_fluid, rho_fluid, d1, d3, msmpc);
    pt(i) = 2*pi*quadgk(fun, 0, thetamax);
    
    % Time it
    if i == 100
        tme = toc/60*length(f)/i;
        fprintf('Estimated time of arrival: %f min\n', tme)
    end
    
    if debug
        debugplots(theta, freq, Phi, T, Ht, E); %#ok<UNRCH>
    end
end

%% Convolve the pulse and the impulse response of the observation point
yt = conv(y, fft(pt, nfft));
tt = (0:length(yt)-1)/fs;

%% Finnished
if saveresults
    dtstr = datestr(now, 'dd-mm-yyyy_HH-MM-SS'); %#ok<UNRCH>
    fprintf('Finnished: %s\n', dtstr)
    outfile = sprintf('wholeshebang_%s_thickness%d.mat', dtstr, d);
    fprintf('Saved to %s\n', outfile)
    save(outfile);
end