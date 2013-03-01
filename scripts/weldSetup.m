%% Do the whole she bang!
debug = false;
saveresults = true;

%% Samplings stuff
fs = 2e6;
nfft = 2^15;
thetamax = pi/2;
f = (0:nfft-1)*fs/nfft;

%% Transducer specs
% Transducer radius
aTx = 25e-3/2;
aRx = 2.5e-3;
% Distance to pipe wall
d1 = 41e-3;
d3 = 41e-3;
% Transducer separation
deltar = 40e-3;

%% Material parameters
rho_fluid = 1000;
v_fluid = 1500;
v_layer = 5850;
fluid1 = struct('v', v_fluid, 'density', rho_fluid);
fluid3 = fluid1;
layer = struct('v', v_layer, 'density', 7850, 'vShear', 3218);

% Steel thickness
d = 8.8e-3;
fres = 0.5*v_layer/d;

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, d);

%% Excitation pulse
tend = 10e-6;
t = (0:1/fs:tend)';

% Start and stop frequencies
f0 = 200e3;
f1 = 800e3;

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
    
    fun = @(xx) shiftedTransducer(xx, freq, aRx, aTx,...
        v_fluid, rho_fluid, d1, d3, model, 0, deltar);
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