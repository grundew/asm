function [pt, f] = computeAsmIntegral(func, asmParams, varargin)
% computeAsmIntegral(func, 'param1', value1, 'param2', value2, ...)
%
% 
% Example:
% p = generateAsmConfig('water', 'steel');
% p.f = linspace(100e3, 500e3, 500);
% [pt, f] = computeAsmIntegral(@integrandFluidSolidFluid_planepiston, p,...
%                               'MaxIntervalCount', 1000);
% figure
% plot(f, abs(pt))

%% Pars default arguments
if nargin < 2
    fprintf('Must have at least to input arguments.\n');
    return;
end

if ~isa(func, 'function_handle')
    fprintf('First input must be a function handle.\n');
    return;
end

if ~asmParams.debug
    pdir = pwd();
    cd(fileparts(mfilename('fullpath')));
    gitStatus = git('status');
    
    isunmodified = ~isempty(regexp(gitStatus, 'modified', 'ONCE'));
    isnotup2date = ~isempty(regexp(gitStatus, 'up-to-date', 'ONCE'));
    if  isunmodified || isnotup2date
        warning('The ASM repos is modified or not up-to-date!')
    else
        asmParams.gitInfo = getGitInfo();
    end
    cd(pdir)
end

%% Integrate over all angles for the point on the axis
[pt, f] = func(asmParams, varargin{:});


if asmParams.savemat
    %% Save results
    dtstr = datestr(now, 'dd_mm_yyyy_HHMMSS');
    fprintf('Finnished: %s\n', dtstr);
    outfilename = generateFilenameString(asmParams, dtstr);
    fprintf('Saved to %s\n', outfilename);
    save(outfilename, 'asmParams', 'pt', 'f');
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