function [V, f] = computeAsmIntegral(func, asmParams, varargin)
% computeAsmIntegral - Computes the integral over the function func.
% 
% [V, f] = computeAsmIntegral(func, asmParams,...
%                              'param1', value1, 'param2', value2, ...)
%
%
% Input:
% func      - Function handle to integral* function
% asmParams - Struct that describes the model, see parseAsmInput
% varargin  - Additional parameter values to quadgk
%
%
% Output:
% V - Frequency response
% f - Frequencies (Hz)
%
% 
% Example:
% p = generateAsmConfig('water', 'steel');
% p.f = linspace(100e3, 500e3, 500);
% [pt, f] = computeAsmIntegral(@integralFluidSolidFluid, p,...
%                               'MaxIntervalCount', 1000);
% figure
% plot(f, abs(pt))


%% Parse default arguments
if nargin < 2
    fprintf('\n Usage: [pt, f] = computeAsmIntegral(func, asmParams)\n\n');
    return;
end

if ~isa(func, 'function_handle')
    fprintf('\n computeAsmIntegral - First input must be a function handle.\n\n');
    return;
end


%% This requires git to be installed
if ~asmParams.debug
    pdir = pwd();
    cd(fileparts(mfilename('fullpath')));
    
    gitStatus = git('status');
    
    isunmodified = ~isempty(regexp(gitStatus, 'modified', 'ONCE'));
    isnotup2date = ~isempty(regexp(gitStatus, 'up-to-date', 'ONCE'));
    if  isunmodified || isnotup2date
        warning('ASM:NOTUPTODATE',...
            'The ASM repos is modified or not up-to-date!')
    elseif strcmp(gitStatus, 'git is not included in the path')
        warning('ASM:NOGIT',...
            'Git not found in path, version information is not included in saved data.')
    else
        asmParams.gitInfo = getGitInfo();
    end
    cd(pdir)
end


%% Integrate over all angles for the point on the axis
[V, f] = func(asmParams, varargin{:});


%% Save results if savemat is true
if asmParams.savemat
    dtstr = datestr(now, 'dd_mm_yyyy_HHMMSS');
    fprintf('Finnished: %s\n', dtstr);
    outfilename = generateFilenameString(asmParams, dtstr);
    fprintf('Saved to %s\n', outfilename);
    save(outfilename, 'asmParams', 'V', 'f');
end


end


function outfilename = generateFilenameString(parameters, dtestr)
% generateFilenameString - Generates a file name based on the asm
%                          parameter struct
%
%
% outfilename = generateFilenameString(parameters, dtestr)
%
% Input:
% parameters - ASM parameter struct
% dtestr     - Current date in string format
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