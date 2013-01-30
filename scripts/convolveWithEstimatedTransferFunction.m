%% Load the transfer function
trans = load('transferfunction.mat');

%% Samplings stuff
fs = 2e6;
nfft = 2^15;
thetamax = pi/2;
f = (0:nfft-1)*fs/nfft;

%% Resample transfer function to current sampling frequency
ht = resample(trans.ht, fs, trans.fs);

%% Transducer specs
aTx = 9e-3;
aRx = 3e-3;
d1 = 10e-2;
d3 = 10e-2;

%% Material parameters
rho_fluid = 1.5;
v_fluid = 350;
% rho_fluid = 1000;
% v_fluid = 1500;

v_layer = 5850;
damping = 1;
fluid1 = struct('v', v_fluid, 'density', rho_fluid);
fluid3 = fluid1;
layer = struct('v', v_layer, 'density', 7850, 'vShear', 3218);
d = 9.8e-3;
fres = 0.5*v_layer/d;

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, d);
thresh = 1e-10;

%% Damping coefficients
alphaLambda = 10.^(0.08./20);

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
    
    fun = @(xx) orofiniIntegrand(xx, freq, aRx, aTx,...
        v_fluid, rho_fluid, d1, d3, model, alphaLambda);
    pt(i) = 2*pi*quadgk(fun, 0, thetamax);
    
    % Time it
    if i == 100
        tme = toc/60*length(f)/i;
        fprintf('Estimated time of arrival: %f min\n', tme)
    end
    
    % if debug
    %    debugplots(theta, freq, Phi, T, Ht, E); %#ok<UNRCH>
    % end
end

%% Convolve the pulse and the impulse response of the observation point
yt = conv(ht, fft(pt, nfft));
tt = (0:length(yt)-1)/fs;