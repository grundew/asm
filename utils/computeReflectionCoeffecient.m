function [X, f, theta] = computeReflectionCoeffecient(ntheta, params)
% [X, f, theta] = computeReflectionCoeffecient(ntheta, params)
%


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


%% Samplings stuff
f = params.f;
nf = length(f);
thetamax = params.thetamax;
thetamin = 0;
theta = linspace(thetamin, thetamax, ntheta);
q = sin(theta);

X = zeros(ntheta, nf);
for i = 1:nf
    %% Plate response, angular
    % Multiply with wave length and convert from dB to linear
    % Loss parameter
    if al_dB ~= 0
        % log(10)/10 = 0.2303
        alphaL = al_dB*0.2303*f(i)/c_F;
    else
        alphaL = 0;
    end
    
    %% Reflection/Transmission coefficient
    if refl
        Plate = reflectionCoefficientAnalytical(f(i), q,...
            thick, rho_F, rho_S, cp, cs, c_F, alphaL);
    else
        [~, Plate] = reflectionCoefficientAnalytical(f(i), q,...
            thick, rho_F, rho_S, cp, cs, c_F, alphaL);
    end
    
    X(:, i) = Plate;
end

end