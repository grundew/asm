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
zTx = params.distanceTx;
al_dB = params.alphaLambda_dB;
gamma = params.gamma;
reflection = params.reflection;
perfReflection = params.perfectReflection;
prgbar = params.progressbar;


%% Samplings stuff
f = params.f;
thetamax = pi/2;
thetamin = -pi/2 + gamma;


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
    
    
    fun = @(theta_z ) integrand(...
        theta_z, f(i), aTx, aRx, zTx, c_F, rho_F,...
    reflection, perfReflection, al_dB, rho_S, cp, cs, thick, gamma);

    X(i) = quadgk(fun, thetamin, thetamax, varargin{:});
    
end


end


function I = integrand(theta_z, f,...
    aRx, aTx, z1, c_F, rho_F,...
    reflection, perfectReflection,...
    al_dB, rho_S, c_Lr, c_Sr, thick, gamma)


%% Angular frequency and total length of wave vector
w = 2*pi*f;
k = w./c_F;


%% Transmitter spatial spectrum
xx = k*sin(theta_z)*aTx;
W = besselj(1, xx)./xx;
W(xx==0) = 0.5;
PhiTx = 2*pi*aTx^2*W;


%% Receiver spatial spectrum
% See confluence page on Angular Spectrum Method for details on the
% relations between the angles
theta_rx = 2*gamma - theta_z;
q_rx = sin(theta_rx);
kr_rx = k*q_rx;
xx = kr_rx*aRx;
W = besselj(1, xx)./xx;
W(xx==0) = 0.5;
PhiRx = 2*pi*aRx^2*W;


%% Plate response, angular
% See confluence page on Angular Spectrum Method for details on the
% relations between the angles
theta_plate = theta_z - gamma;


%% Loss in plate
% Multiply with wave length and convert from dB to linear
% Loss parameter
if al_dB ~= 0
    % log(10)/10 = 0.2303
    alphaL = al_dB*0.2303*f/c_F;
else
    alphaL = 0;
end


%% Reflection/Transmission coefficient
if perfectReflection

    % If perfect reflection the reflection coefficient is just 1
    Plate = ones(size(theta_z));
    
else

    if reflection
        
        % Reflection coefficient is used
        Plate = reflectionCoefficientAnalytical(f, sin(theta_plate),...
            thick, rho_F, rho_S, c_Lr, c_Sr, c_F, alphaL);

    else
        
        % Transmission coefficient is used
        [~, Plate] = reflectionCoefficientAnalytical(f, sin(theta_plate),...
            thick, rho_F, rho_S, c_Lr, c_Sr, c_F, alphaL);
        
    end
    
end


%% Phase shift from transmitter to plate and from plate to receiver
z = z1 + z1*cos(2*gamma);
x = -z1*sin(2*gamma);
Phase = exp(1i*k*(cos(theta_z)*z + sin(theta_z)*x));

I = k*cos(theta_z).*cos(theta_z + 2*gamma).*Plate.*PhiRx.*PhiTx.*Phase;

if any(isnan(I))
    fprintf('NaN value detected at frequency %f and angle %f\n', f, theta_z(isnan(I)));
end

end