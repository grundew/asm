function naturalGas_steel_naturalGas_txAngle2D(alpha, dist)
%% Do the whole she bang!
saveresults = true;

%% Samplings stuff
% fs = 2.5e6;
% nfft = 2^15;
% ntheta = 2^12;
% thetamax = 0.8;
% thetamax = pi/2;
% f = (0:nfft-1)*fs/nfft;
f = 400e3:0.5e3:1.2e6;

%% Transducer specs
% aTx = 1.88e-2;
% aRx = 1.88e-2;
% dist = 7.5e-2;
% d = 0.945e-2;

aTx = 6e-3;
aRx = 6e-3;
d = 25e-3;

%% Material parameters
% rho_fluid = 1000;
% v_fluid = 1500;
% v_layer = 6420;
% rho_layer = 2700;
% vs_layer = 3040;

rho_fluid = 150;
v_fluid = 400;
v_layer = 5900;
rho_layer = 7850;
vs_layer = 3210;
fluid1 = struct('v', v_fluid, 'density', rho_fluid);
fluid3 = fluid1;
layer = struct('v', v_layer, 'density', rho_layer, 'vShear', vs_layer);

fres = 0.5*v_layer/d;

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, d);


%% Integrate over all angles for the point on the axis
nf = length(f);
pt = zeros(nf, 1);
thetazmin = -alpha;
thetazmax = pi/2;

for i = 1:nf
    % Time it
    if i == 1
        fprintf('Started: %s\n', datestr(now, 'dd-mm-yyyy_HH-MM-SS'));
        tic
    end
    
    freq = f(i);

    fun = @(theta_z, theta_x) integrandFluidSolidFluid_2Dangle(theta_z,...
        freq, aRx, aTx,...
        v_fluid, rho_fluid, dist,...
        model, alpha);
    pt(i) = quadgk(fun, thetazmin, thetazmax);
    
    % Time it
    if i == 300
        tme = toc/60*length(f)/i;
        fprintf('Estimated time of arrival: %f min\n', tme)
    end
    
end

%% Finnished
if saveresults
    dtstr = datestr(now, 'dd-mm-yyyy_HH-MM-SS');
    fprintf('Finnished: %s\n', dtstr)
    outfile = sprintf('asm_alpha_%d_%s.mat', alpha*180/pi, dtstr);
    fprintf('Saved to %s\n', outfile)
    save(outfile);
end

end
