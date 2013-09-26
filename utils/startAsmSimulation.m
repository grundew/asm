function startAsmSimulation(varargin)
% startAsmSimulation('param1', value1, 'param2', value2, ...)
% 
%%
% Valid parameters (all of them have default values)
%
% Solid properties:
% 'thickness'
% 'cp'
% 'cs'
% 'rho_solid'
% 'alphaLambda_dB' - Damping in the solid.
%
% Fluid properties:
% 'cf'
% 'rho_fluid'
%
% Geometric setup:
% 'aTx'
% 'aRx'
% 'distanceTx'
% 'distanceRx'
% 'alpha' - Misalignment angle of transducer (not implemented)
%
% Sampling stuff:
% 'fs'
% 'nfft'
% 'thetamax'
%
% Admin stuff:
% 'filenamevars' - Cell array of parameters with value included in the filename

%% Parse the input
params = parseWNIInput(varargin{:});

fluid1 = struct('v', params.cf, 'density', params.rho_fluid);
fluid3 = fluid1;
layer = struct('v', params.cp, 'density', params.rho_solid, 'vShear', params.cs);

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, params.thickness);

%% Samplings stuff
fs = params.fs;
nfft = params.nfft;
thetamax = params.thetamax;
f = (0:nfft-1)*fs/nfft;

%% Excitation pulse
tend = 50e-6;
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
y = conj(hilbert(y));

% Pad with zeros
t = (0:length(y)-1)/fs'; %#ok<NASGU>
Y = ifft(y, nfft); %#ok<NASGU>

%% Unpack parameters
aRx = params.aRx;
aTx = params.aTx;
v_fluid = params.cf;
rho_fluid = params.rho_fluid;
d1 = params.distanceTx;
d3 = params.distanceRx;
alphaLambda_dB = params.alphaLambda_dB;

%% Integrate over all angles for the point on the axis
tic
nf = length(f);
pt = zeros(nf, 1);
for i = 1:nf
    % Time it
    if i == 1
        fprintf('Started: %s\n', datestr(now, 'dd-mm-yyyy_HHMMSS'));
        tic
    end
    
    freq = f(i);
    fun = @(xx) integrandFluidSolidFluidTransmission_withLoss(xx, freq, aRx, aTx,...
       v_fluid, rho_fluid, d1, d3, model, alphaLambda_dB);
    pt(i) = 2*pi*quadgk(fun, 0, thetamax);

    
    % Time it
    if i == 300
        tme = toc/60*length(f)/i;
        fprintf('Estimated time of arrival: %f min\n', tme)
    end

end

%% Convolve the excitation pulse with the system response
yt = conv(y, fft(pt, nfft));
tt = (0:length(yt)-1)/fs; %#ok<NASGU>

%% Save results
dtstr = datestr(now, 'dd_mm_yyyy_HHMMSS');
fprintf('Finnished: %s\n', dtstr);
outfilename = generateFilenameString(params, dtstr);
fprintf('Saved to %s\n', outfilename);
save(outfilename, 'params', 'pt', 't', 'y', 'tt', 'yt');
end

function outfilename = generateFilenameString(parameters, dtestr)
prefix = 'asm';
fnvars = parameters.filenamevars;
c = cellfun(@(x) sprintf('%s-%d', x, parameters.(x)) , fnvars, 'uni', 0);
paramstr = strjoin(c, '-');
outfilename = sprintf('%s-%s-%s.mat', prefix, paramstr, dtestr);
end
