function I = computeAsmIntegrand(f, theta, varargin)
% I = computeAsmIntegrand(f, theta, varargin)
% 
%
% Valid parameters (all of them have default values)
%
% Solid properties:
% 'thickness'
% 'cp'
% 'cs'
% 'rho_solid'
% 'alphaLambda_dB' - Damping in the solid.
%
% Fluid properties:
% 'cf'
% 'rho_fluid'
%
% Geometric setup:
% 'aTx'
% 'aRx'
% 'distanceTx'
% 'distanceRx'
% 'alpha' - Misalignment angle of transducer (not implemented)
%

%% Parse the input
[params, p] = parseAsmInput(varargin{:});

invalidparams = {'fs', 'thetamax', 'nfft', 'filenamevars'};
invid = cellfun(@(x) ~any(strcmp(x, p.UsingDefaults)), invalidparams);
if any(invid)
    error('InputError:InvalidParameter', 'Invalid parameter specified: %s', invalidparams{~invid});
end

fluid1 = struct('v', params.cf, 'density', params.rho_fluid);
fluid3 = fluid1;
layer = struct('v', params.cp, 'density', params.rho_solid, 'vShear', params.cs);

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, params.thickness);

%% Unpack parameters
aRx = params.aRx;
aTx = params.aTx;
v_fluid = params.cf;
rho_fluid = params.rho_fluid;
d1 = params.distanceTx;
d3 = params.distanceRx;
alphaLambda_dB = params.alphaLambda_dB;

%% Integrate over all angles for the point on the axis
nf = length(f);
ntheta = length(theta);
I = zeros(nf, ntheta);
for i = 1:nf
    freq = f(i);
    I(i, :) = integrandFluidSolidFluidTransmission_withLoss(theta, freq, aRx, aTx,...
       v_fluid, rho_fluid, d1, d3, model, alphaLambda_dB);
end

end