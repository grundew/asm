function [pt, f] = air_steel_air_displacement_pointrx(x0, savemat)

%% Samplings stuff
fs = 5e6;
nfft = 2^17;
thetamin = 0;
thetamax = 0.8;
f = (0:nfft-1)*fs/nfft;

%% Transducer specs
aTx = 9e-3;
d_z = 200e-3;
thickness = 11.85e-3;

%% Material parameters
rho_fluid = 1.5;
v_fluid = 342;
v_layer = 5970;
rho_layer = 7850;
vs_layer = 3200;

fluid1 = struct('v', v_fluid, 'density', rho_fluid);
fluid3 = fluid1;
layer = struct('v', v_layer, 'density', rho_layer, 'vShear', vs_layer);

fres = 0.5*v_layer/thickness;

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, thickness);
reflection = false;

%% Integrate over all angles for the point on the axis
nf = length(f);
pt = zeros(nf, 1);

for i = 1:nf/2+1
    % Time it
    if i == 1
        fprintf('Started: %s\n', datestr(now, 'dd-mm-yyyy_HH-MM-SS'));
        tic
    end
    
    freq = f(i);
    fun = @(theta) integrandFluidSolidFluid_pointrx(theta, freq, aTx,...
        v_fluid, rho_fluid, d_z, x0, model, reflection);
    pt(i) = k*quadgk(fun, thetamin, thetamax);
    
    % Time it
    if i == 300
        tme = toc/60*length(f)/i;
        fprintf('Estimated time of arrival: %f min\n', tme)
    end
    
end
k = 2*pi*f'/v_fluid;
pt = 4*pi/rho_fluid/v_fluid/aTx*k.*pt;

%% Finnished
dtstr = datestr(now, 'dd_mm_yyyy_HHMMSS');
fprintf('Finnished: %s\n', dtstr)
if exist('savemat', 'var') && savemat
    outfile = sprintf('asm_tt_displaceRx_%d_%s.mat', x0*1e3, dtstr);
    fprintf('Saved to %s\n', outfile)
    save(outfile);
end

end