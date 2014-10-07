function [pt, f, params] = computeAsmIntegral(func, varargin)
% startAsmSimulation(func, 'param1', value1, 'param2', value2, ...)
% 
% Input:
% func - Cell array. First element must be a function handle and the rest
%        of the elements are the non-default arguments to the function. 
%        The function must take the following arguments (but might ignore them):
%          1. Angle of incidence for the plane wave (radians)
%          2. Frequency (Hz)
%          3. Radius of receiver (m)
%          4. Radius of transmitter (m)
%          5. Speed of sound in fluid (m/s)
%          6. Distance from transmitter to target (m)
%          7. Distance from receiver to target (m)
%          8. Model object/struct
%          9. Loss in steel, $\alpha \lambda$ (dB)
%          10. Reflection, true if reflection, false if transmission
% 
% Valid parameters (all of them have default values)
%
% Solid properties:
% 'thickness'
% 'cp'
% 'cs'
% 'rho_solid'
% 'alphaLambda_dB' - Damping in the solid.
% 'reflection' - (Boolean) Reflection coefficient is used if true
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
% 'displaceRx' - Displacement of receiver perpendicular to the acuostical
%                axis (m)
%
% Sampling stuff:
% 'f'
% 'thetamin' - Integration lower limit
% 'thetamax' - Integration limit
%
% Admin stuff:
% 'filenamevars' - Cell array of parameters with value included in the filename

%% Pare default arguments
if isa(func, 'function_handle')
    func_args = {};
    func_handle = func;
elseif iscell(func) && isa(func{1}, 'function_handle')
    if length(func) > 1
        func_args = func{2:end};
    else
        func_args = {};
    end
    func_handle = func{1};
else
    
end

%% Parse the parameters
params = parseAsmInput(varargin{:});


%% Integrate over all angles for the point on the axis
pt = func_handle(params, func_args{:});


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
idf = strcmp(fnvars, 'f');

if any(idf)
    fnvars = {fnvars{~idf}, 'fmin', 'fmax'};
    parameters.('fmin') = min(parameters.f);
    parameters.('fmax') = max(parameters.f);
    df = diff(parameters.f);
    if length(unique(df)) == 1
        parameters.df = unique(df);
        fnvars{end+1} = 'df';
    end
end
    
if isempty(fnvars)
    outfilename = sprintf('%s_%s.mat', prefix, dtestr);
else
    c = cellfun(@(x) sprintf('%s_%d', x, parameters.(x)) , fnvars, 'uni', 0);
    paramstr = strjoin(c, '_');
    outfilename = sprintf('%s_%s_%s.mat', prefix, paramstr, dtestr);
end

end