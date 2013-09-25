%% Do the whole she bang!
debug = false;
saveresults = true;

%% Samplings stuff
fs = 2.5e6;
nfft = 2^15;
ntheta = 2^12;
thetamax = 0.8;
% thetamax = pi/2;

f = (0:nfft-1)*fs/nfft;

%% Transducer specs
aTx = 6e-3;
aRx = 6e-3;
dist = 0.08;

%% Material parameters
rho_fluid = 120;
v_fluid = 430;

v_layer = 5850;
fluid1 = struct('v', v_fluid, 'density', rho_fluid);
fluid3 = fluid1;
layer = struct('v', v_layer, 'density', 7850, 'vShear', 3218);
d = 25e-3;
fres = 0.5*v_layer/d;

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, d);
thresh = 1e-10;

% Misalignment angle
% alpha = 1.5/180*pi;
alpha = 0;

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
theta = linspace(0, thetamax, ntheta);
dtheta = theta(2) - theta(1);

for i = 1:nf
    % Time it
    if i == 1
        fprintf('Started: %s\n', datestr(now, 'dd-mm-yyyy_HH-MM-SS'));
        tic
    end
    
    freq = f(i);

    fun = @(xx) integrandFluidSolidFluid_withAngle(xx, freq, aRx, aTx,...
       v_fluid, rho_fluid, dist, model, alpha);
    pt(i) = 2*pi*quadgk(fun, 0, thetamax);
    % Ht = integrandFluidSolidFluid_withAngle(theta, freq, aRx, aTx,...
    %     v_fluid, rho_fluid, dist, model, alpha);
    % pt(i) = 2*pi*2*dtheta*(0.5*(Ht(1) + Ht(end)) + sum(Ht(2:end-1)));
    
    % Time it
    if i == 300
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