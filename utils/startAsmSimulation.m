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
params = parseInput(varargin{:});

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

function result = parseInput(varargin)
%% Default values
% Solid
def_thickness = 10e-3; % m
def_cs = 3158; % m/s
def_cp = 5850; % m/s
def_rho_solid = 7850; % m/s
def_alphaLambda = 0.008; % dB
% Fluid
def_cf = 342.21; % m/s
def_rho_fluid = 1.5; % kg/m^3
% Geometric setup
def_aTx = 9e-3; % Radius (m)
def_aRx = 3e-3; % Radius (m)
def_distanceTx = 40e-3; % m
def_distanceRx = 40e-3; % m
% Sampling stuff
def_fs = 2e6;
def_nfft = 2^15;
def_thetamax = 0.8;

%% Parse input
p = inputParser();

% Validator, true if x is a scalar greater than zero (gtz)
validatorfun = @(x) validateattributes(x,...
    {'numeric'}, {'scalar', 'real', 'nonempty', '>', 0});
validatorfunzero = @(x) validateattributes(x,...
    {'numeric'}, {'scalar', 'real', 'nonempty', '>=', 0});

% Solid
p.addParamValue('thickness', def_thickness, validatorfun)
p.addParamValue('cp', def_cp, validatorfun);
p.addParamValue('cs', def_cs, validatorfun);
p.addParamValue('rho_solid', def_rho_solid, validatorfun);
p.addParamValue('alphaLambda_dB', def_alphaLambda, validatorfunzero);
% Fluid
p.addParamValue('cf', def_cf, validatorfun);
p.addParamValue('rho_fluid', def_rho_fluid, validatorfun);
% Geometric setup
p.addParamValue('aTx', def_aTx, validatorfun);
p.addParamValue('aRx', def_aRx, validatorfun);
p.addParamValue('distanceTx', def_distanceTx, validatorfun);
p.addParamValue('distanceRx', def_distanceRx, validatorfun);
% Sampling stuff
p.addParamValue('fs', def_fs, validatorfun);
p.addParamValue('nfft', def_nfft, validatorfun);
p.addParamValue('thetamax', def_thetamax, validatorfun);
% Admin stuff
p.addParamValue('filenamevars', {}, @iscell);
% Parse
p.parse(varargin{:});
result = p.Results;
end

function outfilename = generateFilenameString(parameters, dtestr)
prefix = 'asm';
fnvars = parameters.filenamevars;
c = cellfun(@(x) sprintf('%s-%d', x, parameters.(x)) , fnvars, 'uni', 0);
paramstr = strjoin(c, '-');
outfilename = sprintf('%s-%s-%s.mat', prefix, paramstr, dtestr);
end
