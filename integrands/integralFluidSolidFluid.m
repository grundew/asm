function [X, f] = integralFluidSolidFluid(params, varargin)
% integrandFluidSolidFluidPlanePistion - computes the integral for an axial
% symmetric model with either a point, plane or focused source.
%
% [X, f] = integrandFluidSolidFluidPlanePiston(params, varargin)
%
% Output:
% X - Computed integral at frequencies, f
% f - Frequencies
%
% References:
% 1. Orofino, 1992. http://dx.doi.org/10.1121/1.405408


%% Samplings stuff
f = params.f;
thetamax = params.thetamax;
thetamin = 0;

%% Unpack parameters
aRx = params.aRx;
aTx = params.aTx;
c_F = params.cf;
rho_F = params.rho_fluid;
rho_S = params.rho_solid;
cp = params.cp;
cs = params.cs;
thick = params.thickness;
zTx = params.distanceTx;
zRx = params.distanceRx;
al_dB = params.alphaLambda_dB;
x0 = params.displaceRx;
reflection = params.reflection;
perfReflection = params.perfectReflection;
prgbar = params.progressbar;
focusRx = params.focusRx;
focusTx = params.focusTx;


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
    
    
    fun = @(theta_z ) integrandFluidSolidFluidAxialSymmetric(...
        theta_z, f(i), k, zTx, zRx, aTx, aRx, c_F, rho_F, x0, focusRx, focusTx,...
    reflection, perfReflection, al_dB, rho_S, cp, cs, thick);

    X(i) = quadgk(fun, thetamin, thetamax, varargin{:});
    
end

end