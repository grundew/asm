function [X, f, theta, p] = computeReflectionCoeffecient(ntheta, varargin)
% X = computeReflectionCoeffecient(ntheta, varargin)
%

%% Parse the parameters
p = parseAsmInput(varargin{:});

fluid1 = struct('v', p.cf, 'density', p.rho_fluid);
fluid3 = fluid1;
layer = struct('v', p.cp, 'density', p.rho_solid, 'vShear', p.cs);

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, p.thickness);

%% Samplings stuff
f = p.f;
nf = length(f);
thetamax = p.thetamax;
thetamin = 0;
theta = linspace(thetamin, thetamax, ntheta);

X = zeros(ntheta, nf);
for i = 1:nf
    %% Reflection/Transmission coefficient
    if p.reflection
        Plate = analyticRTFast(f(i), theta, model);
    else
        Plate = transmissionCoefficientAnalytical(f(i), sin(theta), model);
    end
    
    X(:, i) = Plate;
end

end