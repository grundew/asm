function [X, f] = integralFluidSolidFluid_pointrx(params, varargin)
% [X, f] = integrandFluidSolidFluid_pointrx(params varargin)
%
% Output:
% X - Computed integral at frequencies f
% f - Frequencies
%


%% Samplings stuff
f = params.f;
thetamax = params.thetamax;
thetamin = 0;

%% Unpack parameters
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
prgbar = params.progressbar;


%% Check sanity of parameters and give warnings or errors
% TODO: Implement sanity checks. Throw errors.
assert(al_dB >=0, 'asm:paramerror', 'The damping must be a positive number.');


% ASCII Progress bar
if prgbar
    nnn  = progress('init', 'Wait');
    tic       % only if not already called
    t0 = toc; % or toc if tic has already been called
    tm = t0;
end

X = zeros(size(f));
nf = length(f);
for i = 1:nf
    % ASCII Progress bar
    if prgbar
        tt = ceil((toc-t0)*(nf-i)/i);
        nn = datenum(0, 0, 0, 0, 0, tt);
        nnvec = datevec(nn);
        prgstr = sprintf('ETA: %.0f HOURS %.0f MIN %.0f SEC',...
            nnvec(end-2), nnvec(end-1), ceil(nnvec(end)));
        progress(i/nf, sprintf('%i/%i (%s)', i, nf, prgstr));
    end
    
    w = 2*pi*f(i);
    k = w/c_F;
    fun = @(theta_z) integrandFluidSolidFluid_pointrx(...
        theta_z, f(i), k, d1, d3, aTx, c_F,...
        al_dB, refl, rho_F, rho_S, cp, cs, thick, x0);
    X(i) = quadgk(fun, thetamin, thetamax, varargin{:});
    
end

end