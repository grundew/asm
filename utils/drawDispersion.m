function drawDispersion(f, theta, varargin)
% drawDispersion(f, theta, 'param1', value1, 'param2', value2, ...)
%
% Input:
% f - Frequency (Hz)
% theta - Angle (rad)
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

%% Parse the input
params = parseAsmInput(varargin{:});

% invalidparams = {'fs', 'thetamax', 'distanceTx', 'distanceRx',...
%     'aTx', 'aRx', 'filenamevars', 'nfft'};
% invid = cellfun(@(x) ~any(strcmp(x, p.UsingDefaults)), invalidparams);
% if any(invid)
%     errorstr = sprintf('Invalid parameter specified: %s', invalidparams{~invid});
%     error('InputError:InvalidParameter', errorstr);
% end

fluid1 = struct('v', params.cf, 'density', params.rho_fluid);
fluid3 = fluid1;
layer = struct('v', params.cp, 'density', params.rho_solid, 'vShear', params.cs);

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, params.thickness);

fres = params.cp/params.thickness/2;

% Calculate the reflection coefficient for all f and theta
R = zeros(length(f), length(theta));
for i = 1:length(f)
    R(i, :) = reflectionTransmissionCoffecientAnalytical(f(i), theta, model);
end

%% Solid layer immersed in a fluid
plotdisp(f/fres, theta*180/pi, db(abs(R))')
end

function plotdisp(x, y, R)
%% Plot
fig = figure;
imagesc(x, y, R)
axis xy
xlabel('f/f_1')
ylabel('\theta (deg)')
colormap(gray)
end