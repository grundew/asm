function [pt, f, params] = startAsmSimulation(varargin)
% startAsmSimulation('param1', value1, 'param2', value2, ...)
% 
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
% 'f'
% 'thetamax'
%
% Admin stuff:
% 'filenamevars' - Cell array of parameters with value included in the filename

%% Parse the input
params = parseAsmInput(varargin{:});

fluid1 = struct('v', params.cf, 'density', params.rho_fluid);
fluid3 = fluid1;
layer = struct('v', params.cp, 'density', params.rho_solid, 'vShear', params.cs);

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, params.thickness);

%% Samplings stuff
f = params.f;
thetamax = params.thetamax;

%% Unpack parameters
aRx = params.aRx;
aTx = params.aTx;
v_fluid = params.cf;
rho_fluid = params.rho_fluid;
d1 = params.distanceTx;
d3 = params.distanceRx;
alphaLambda_dB = params.alphaLambda_dB;
fres = 0.5*params.cp/params.thickness; %#ok<*NASGU>

%% Integrate over all angles for the point on the axis
tic
nf = numel(f);
pt = zeros(size(f));
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

if params.savemat
    %% Save results
    dtstr = datestr(now, 'dd_mm_yyyy_HHMMSS');
    fprintf('Finnished: %s\n', dtstr);
    outfilename = generateFilenameString(params, dtstr);
    fprintf('Saved to %s\n', outfilename);
    save(outfilename, 'params', 'pt', 'f', 'fres');
end
end

function outfilename = generateFilenameString(parameters, dtestr)
prefix = 'asm';
fnvars = parameters.filenamevars;
c = cellfun(@(x) sprintf('%s_%d', x, parameters.(x)) , fnvars, 'uni', 0);
paramstr = strjoin(c, '_');
outfilename = sprintf('%s_%s_%s.mat', prefix, paramstr, dtestr);
end