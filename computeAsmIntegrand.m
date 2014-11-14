function [X, theta] = computeAsmIntegrand(func, ntheta, params)
% [X, theta] = computeAsmIntegrand(func, ntheta, asmParams)
%
% 
% Example:
% p = generateAsmConfig('water', 'steel');
% p.f = linspace(100e3, 500e3, 500);
% ntheta = 500;
% [X, theta] = computeAsmIntegrand(...
%                 @integrandFluidSolidFluid_planepiston, ntheta, p);
%
% figure
% imagesc(p.f, theta, abs(X))

%% Pars default arguments
if nargin < 3
    fprintf('Must have at least to input arguments.\n');
    return;
end

if ~isa(func, 'function_handle')
    fprintf('First input must be a function handle.\n');
    return;
end

if ~params.debug
    pdir = pwd();
    cd(fileparts(mfilename('fullpath')));
    gitStatus = git('status');
    
    isunmodified = ~isempty(regexp(gitStatus, 'modified', 'ONCE'));
    isnotup2date = ~isempty(regexp(gitStatus, 'up-to-date', 'ONCE'));
    if  isunmodified || isnotup2date
        warning('The ASM repos is modified or not up-to-date!')
    else
        params.gitInfo = getGitInfo();
    end
    cd(pdir)
end


%% Samplings stuff
f = params.f;
thetamax = params.thetamax;
thetamin = 0;
theta = linspace(thetamin, thetamax, ntheta);


%% Unpack parameters
aRx = params.aRx;
aTx = params.aTx;
c_F = params.cf;
rho_F = params.rho_fluid;
rho_S = params.rho_solid;
cp = params.cp;
cs = params.cs;
thick = params.thickness;
d1 = params.distanceTx;
d3 = params.distanceRx;
al_dB = params.alphaLambda_dB;
fres = 0.5*params.cp/params.thickness; %#ok<*NASGU>
x0 = params.displaceRx;
refl = params.reflection;


%% Check sanity of parameters and give warnings or errors
% TODO: Implement sanity checks. Throw errors.
assert(al_dB >=0, 'asm:paramerror', 'The damping must be a positive number.');


%% Evaluate integrand over all angles for the point on the axis
X = zeros(length(f), ntheta);
for i = 1:length(f)    
    w = 2*pi*f(i);
    k = w/c_F;

    X(i, :) = func(theta, f(i), k, d1, d3, aTx, aRx, c_F,...
        al_dB, refl, rho_F, rho_S, cp, cs, thick, x0);
end

end