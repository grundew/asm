function [pt, f, params] = startAsmWithMisalignment_2D(alpha_plate, varargin)
% startAsmSimulation('param1', value1, 'param2', value2, ...)
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
refl = params.reflection;

%% Integrate over all angles for the point on the axis
nf = numel(f);
pt = zeros(size(f));
thetazmin = -pi/2;
thetazmax = pi/2;

for i = 1:nf
    % Time it
    if i == 1
        fprintf('Started: %s\n', datestr(now, 'dd-mm-yyyy_HH-MM-SS'));
        tstart = tic;
    end
    
    freq = f(i);

    fun = @(theta_z) integrandFluidSolidFluid_2Dangle(theta_z,...
        freq, aRx, aTx,...
        v_fluid, rho_fluid, d1+d3,...
        model, alpha_plate, refl);
    pt(i) = quadgk(fun, thetazmin, thetazmax);
    
    % Time it
    if i == 300
        tme = toc(tstart)*length(f)/i;
        nn = datenum(0, 0, 0, 0, 0, tme);
        nnvec = datevec(nn);
        fprintf('Estimated time of arrival: %.0f HOURS %.0f MIN %.0f SEC\n',...
            nnvec(end-2), nnvec(end-1), ceil(nnvec(end)));
    end    
end


if params.savemat
    %% Save results
    dtstr = datestr(now, 'dd_mm_yyyy_HHMMSS');
    fprintf('Finnished: %s\n', dtstr);
    outfilename = generateFilenameString(alpha_plate, dtstr);
    fprintf('Saved to %s\n', outfilename);
    save(outfilename, 'params', 'pt', 'f', 'fres');
end

end

function outfilename = generateFilenameString(alpha_plate, dtestr)
prefix = 'asm';
paramstr = sprintf('alphaplate_%.2f', alpha_plate*180/pi);
outfilename = sprintf('%s_%s_%s.mat', prefix, paramstr, dtestr);
end
