function water_steel_water_displacement_pointrx(x0)

%% Samplings stuff
fs = 5e6;
nfft = 2^15;
ntheta = 2^12;
thetamax = pi/2;
f = (0:nfft-1)*fs/nfft;
theta = linspace(0, thetamax, ntheta);

%% Transducer specs
aTx = 6e-3;
d = 102e-3;
thickness = 11.85e-3;

%% Material parameters
rho_fluid = 1000;
v_fluid = 1450;
v_layer = 5970;
rho_layer = 7850;
vs_layer = 3210;

fluid1 = struct('v', v_fluid, 'density', rho_fluid);
fluid3 = fluid1;
layer = struct('v', v_layer, 'density', rho_layer, 'vShear', vs_layer);

fres = 0.5*v_layer/d;

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, thickness);
reflection = false;

%% Integrate over all angles for the point on the axis
nf = length(f);
pt = zeros(nf, 1);

for i = 1:nf
    % Time it
    if i == 1
        fprintf('Started: %s\n', datestr(now, 'dd-mm-yyyy_HH-MM-SS'));
        tic
    end
    
    freq = f(i);
    I = integrandFluidSolidFluid_pointrx(theta, freq, aTx,...
        v_fluid, rho_fluid, d, x0, model, reflection);
    pt(i) = trapz(theta, I);
    % pt(i) = quadgk(fun, thetazmin, thetazmax);
    
    % Time it
    if i == 300
        tme = toc/60*length(f)/i;
        fprintf('Estimated time of arrival: %f min\n', tme)
    end
    
end

%% Finnished
dtstr = datestr(now, 'dd_mm_yyyy_HHMMSS');
fprintf('Finnished: %s\n', dtstr)
outfile = sprintf('asm_displaceRx_%d_%s.mat', x0*1e-3, dtstr);
fprintf('Saved to %s\n', outfile)
save(outfile);

end
