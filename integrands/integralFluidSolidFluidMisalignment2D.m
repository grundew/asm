function [X, f] = integralFluidSolidFluidMisalignment2D(params, varargin)
% [X, f] = integralFluidSolidFluidMisalignment2D(params, varargin)


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
fres = 0.5*params.cp/params.thickness; %#ok<*NASGU>
refl = params.reflection;
alpha = params.alpha;
prgbar = params.progressbar;


%% Samplings stuff
f = params.f;
thetamax = params.thetamax;
thetamin = -alpha;


% ASCII Progress bar
if prgbar
    nnn  = progress('init', 'Wait');
    tic       % only if not already called
    t0 = toc; % or toc if tic has already been called
    tm = t0;
end


nf = length(f);
X = zeros(size(f));
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
    fun = @(theta_z) integrand(...
        theta_z, f(i), k, d1, aTx, aRx, c_F,...
        refl, rho_F, rho_S, cp, cs, thick);
    X(i) = quadgk(fun, thetamin, thetamax, varargin{:});
    
end


end


function I = integrand(theta_z, f, k, d1,...
   aTx, aRx, c_F, refl, rho_F, rho_S, c_Lr, c_Sr, thick)
%% Angular frequency and total length of wave vector
w = 2*pi*f;


%% Transmitter spatial spectrum
xx = k*sin(theta_z)*aTx;
W = besselj(1, xx)./xx;
W(xx==0) = 0.5;
PhiTx = 2*pi*aTx^2*W;


%% Receiver spatial spectrum
theta_rx = 2*alpha - theta_z;
q_rx = sin(theta_rx);
kr_rx = k*q_rx;
xx = kr_rx*aRx;
W = besselj(1, xx)./xx;
W(xx==0) = 0.5;
PhiRx = 2*pi*aRx^2*W;


%% Reflection/Transmission coefficient
theta_plate = alpha - theta_z;
q = sin(theta_plate);
if refl
    Plate = reflectionCoefficientAnalytical(f, q,...
        thick, rho_F, rho_S, c_Lr, c_Sr, c_F, alphaL);
else
    [~, Plate] = reflectionCoefficientAnalytical(f, q,...
        thick, rho_F, rho_S, c_Lr, c_Sr, c_F, alphaL);
end


%% Phase shift from transmitter to plate and from plate to receiver
z = d1 + d1*cos(2*alpha);
x = -d1*sin(2*alpha);
Phase = exp(1i*k*(cos(theta_z)*z + sin(theta_z)*x));


%% Final integrand
I = k*cos(theta_z).*cos(theta_z-2*alpha).*Plate.*PhiRx.*PhiTx.*Phase;

if any(isnan(I))
    fprintf('NaN value detected at frequency %f and angle %f\n', f, theta_z(isnan(I)));
end

end
