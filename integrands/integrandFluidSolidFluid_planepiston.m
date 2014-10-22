function [X, f] = integrandFluidSolidFluid_planepiston(params, varargin)
% [X, f] = integrandFluidSolidFluid_planepiston(params varargin)
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
    fun = @(theta_z) integrand(theta_z, f(i), k, d1, d3, aTx, aRx, c_F,...
        al_dB, refl, rho_F, rho_S, cp, cs, thick, x0);
    
    
    X(i) = quadgk(fun, thetamin, thetamax, varargin{:});
    
    % Time it
    if i == 500
        tme = toc(tstart)*length(f)/i;
        nn = datenum(0, 0, 0, 0, 0, tme);
        nnvec = datevec(nn);
        fprintf('Estimated time of arrival: %.0f HOURS %.0f MIN %.0f SEC\n',...
            nnvec(end-2), nnvec(end-1), ceil(nnvec(end)));
    end

end

end

function I = integrand(theta_z, f, k, d1, d3, aTx, aRx, c_F,...
    al_dB, refl, rho_F, rho_S, c_Lr, c_Sr, thick, x0)
%% Compute wave numbers
q = sin(theta_z);
p = sqrt(1-q.^2);

kr = k*q;
kz = k*p;

%% Transmitter spatial sensitivity
% No angle adjustment
xx = kr*aTx;
W = besselj(1, xx)./xx;
W(xx==0) = 0.5;
PhiTx = 2*pi*aTx^2*W;


%% Receiver spatial sensitivity
xx = kr*aRx;
Wrx = besselj(1, xx)./xx;
Wrx(xx==0) = 0.5;
PhiRx = 2*pi*aRx^2*Wrx;


%% Plate response, angular
% Multiply with wave length and convert from dB to linear
% Loss parameter
if al_dB ~= 0
    % log(10)/10 = 0.2303
    alphaL = al_dB*0.2303*f/c_F;
else
    alphaL = 0;
end


%% Reflection/Transmission coefficient
if refl
    Plate = reflectionCoefficientAnalytical(f, q,...
        thick, rho_F, rho_S, c_Lr, c_Sr, c_F, alphaL);
else
    [~, Plate] = reflectionCoefficientAnalytical(f, q,...
        thick, rho_F, rho_S, c_Lr, c_Sr, c_F, alphaL);
end


%% Displacement factor
if x0 > 0
    dispRx = besselj(0, x0*kr);
else
    dispRx = 1;
end


%% Phase shift from transmitter to plate and from plate to receiver
Phase = exp(1i*kz*(d1 + d3));


%% Assemble integrand
I = Plate.*k.*q.*dispRx.*PhiRx.*PhiTx.*Phase.*k.*p.^2;

end