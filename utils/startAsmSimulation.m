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
% 'aTx' - Radius of transmitter
% 'aRx' - Radius of receiver
% 'distanceTx' - Distance from the transmitter to the plate
% 'distanceRx' - Distance from receiver to the plate
% 'alpha' - Misalignment angle of transducer (not implemented)

%
% Sampling stuff:
% 'f'
% 'thetamax' - Integration limit
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
x0 = params.displaceRx;
alpha_plate = params.alpha_plate;

%% Integrate over all angles for the point on the axis
tic
refl = params.reflection;
displaced = x0 ~= 0;
% alpha = 0;

nf = numel(f);
pt = zeros(size(f));

for i = 1:nf
    % Time it
    if i == 1
        fprintf('Started: %s\n', datestr(now, 'dd-mm-yyyy_HHMMSS'));
        tic
    end
    
    freq = f(i);
    w = 2*pi*f(i);
    k = w/v_fluid;
    
    if refl && displaced
        % Reflection and displaced receiver
        % Integrate over k_r from 0 to k
        fun = @(xx) integrandFluidSolidFluidReflection_withLossAndDisplacement(xx, w, aRx, aTx,...
            v_fluid, rho_fluid, d1, d3, model, x0, alphaLambda_dB);
        pt(i) = quadgk(fun, 0, k);
    elseif refl && ~displaced
        % Reflection in pulse-echo mode with same rx and tx
        fun = @(xx) integrandFluidSolidFluidReflection(xx, freq, aRx, aTx,...
            v_fluid, rho_fluid, d1, d3, model);
        pt(i) = 2*pi*quadgk(fun, 0, thetamax);
    elseif ~refl && ~displaced
        % Transmission and rx and tx on same axis
        fun = @(xx) integrandFluidSolidFluidTransmission_withLoss(xx, freq, aRx, aTx,...
            v_fluid, rho_fluid, d1, d3, model, alphaLambda_dB);
        pt(i) = 2*pi*quadgk(fun, 0, thetamax);
    elseif ~refl && displaced
        % Transmission with receiver displaced
        fun = @(xx) integrandFluidSolidFluidTransmission_withLossAndDisplacement(xx, w, aRx, aTx,...
            v_fluid, rho_fluid, d1, d3, model, x0, alphaLambda_dB);
        pt(i) = quadgk(fun, 0, k);
    elseif alpha_plate ~= 0
        fun = @(xx) integrandFluidSolidFluid_withAngle(xx, freq, aRx, aTx,...
            v_fluid, rho_fluid, d1, model, alpha);
        pt(i) = 2*pi*quadgk(fun, 0, thetamax);
    else
        error('Something went haywire');
    end
    
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
if isempty(fnvars)
    outfilename = sprintf('%s_%s.mat', prefix, dtestr);
else
    c = cellfun(@(x) sprintf('%s_%d', x, parameters.(x)) , fnvars, 'uni', 0);
    paramstr = strjoin(c, '_');
    outfilename = sprintf('%s_%s_%s.mat', prefix, paramstr, dtestr);
end

end