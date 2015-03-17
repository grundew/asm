function [X, f] = integralPerfectReflectorPointRx(params, varargin)
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
d1 = params.distanceTx;
d3 = params.distanceRx;
x0 = params.displaceRx;
prgbar = params.progressbar;


%% Check sanity of parameters and give warnings or errors
% TODO: Implement sanity checks. Throw errors.



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
    fun = @(theta_z) integrandPerfectReflectorPointRx(...
        theta_z, w, aTx, c_F, rho_F, d1, d3, x0);
    X(i) = quadgk(fun, thetamin, thetamax, varargin{:});
    
end

end