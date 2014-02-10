function [I, f, theta, params] = calcIntegrand(delta_theta, varargin)
% startAsmSimulation(delta_theta, 'param1', value1, 'param2', value2, ...)
%
% Integrand is calculated for theta = thetamin:delta_theta:theta_max;
% 
% Valid parameters (all of them have default values)
%
% Solid properties:
% 'thickness'
% 'cp'
% 'cs'
% 'rho_solid'
% 'alphaLambda_dB' - Damping in the solid.
% 'alpha' - Misalignment angle of transducer (not implemented)
% 'reflection' - (Boolean) Reflection coefficient is used if true
%
% Fluid properties:
% 'cf'
% 'rho_fluid'
%
% Geometric setup:
% 'aTx' - Radius of transmitter
% 'aRx' - Radius of receiver
% 'distanceTx' - Distance from the transmitter to the plate
% 'distanceRx' - Distance from receiver to the plate
%
% Sampling stuff:
% 'f'
% 'thetamax' - Integration limit
%
% Admin stuff:
% 'filenamevars' - Cell array of parameters with value included in the filename

%% Parse the input
params = parseAsmInput(varargin{:});

fluid1 = struct('v', params.cf, 'density', params.rho_fluid);
fluid3 = fluid1;
layer = struct('v', params.cp, 'density', params.rho_solid, 'vShear', params.cs);

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, params.thickness);

%% Samplings stuff
f = params.f;
thetamax = params.thetamax;

%% Unpack parameters
aRx = params.aRx;
aTx = params.aTx;
v_fluid = params.cf;
rho_fluid = params.rho_fluid;
d1 = params.distanceTx;
d3 = params.distanceRx;
alphaLambda_dB = params.alphaLambda_dB;
fres = 0.5*params.cp/params.thickness; %#ok<*NASGU>
x0 = params.displaceRx;
alpha_plate = params.alpha_plate;
refl = params.reflection;
tx_focus = params.tx_focus;

%% Integrate over all angles for the point on the axis
nf = numel(f);
pt = zeros(size(f));

if abs(alpha_plate) > 0
    thetamin = -thetamax;
else
    thetamin = 0;
end

theta = thetamin:delta_theta:thetamax;
ntheta = length(theta);
I = zeros(nf, ntheta);

for i = 1:nf
    % Time it
    if i == 1
        fprintf('Started: %s\n', datestr(now, 'dd-mm-yyyy_HHMMSS'));
        tstart = tic();
    end
    
    freq = f(i);
    w = 2*pi*f(i);
    k = w/v_fluid;
    
    if ~tx_focus
        % Plane piston with out focus
        I(i, :) = integrandFluidSolidFluid_planepiston(...
            theta, freq, aRx, aTx,...
            v_fluid, d1, d3, model, x0, alphaLambda_dB, refl, alpha_plate);
    elseif tx_focus
        % Plane piston and focused receiver
        I(i, :) = integrandFluidSolidFluid_focusedtx(...
            xx, freq, aRx, aTx,...
            v_fluid, d1, d3,...
            model, x0, alphaLambda_dB, refl, alpha_plate);
    else
        error('Something went haywire');
    end        
    
    % Time it
    if i == 300
        tme = toc(tstart)*length(f)/i;
        nn = datenum(0, 0, 0, 0, 0, tme);
        nnvec = datevec(nn);
        fprintf('Estimated time of arrival: %.0f HOURS %.0f MIN %.0f SEC\n',...
            nnvec(end-2), nnvec(end-1), ceil(nnvec(end)));
    end

end

end