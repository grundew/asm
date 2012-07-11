function drawDispersion()

nt = 100;
thetamax = 0.5;
theta = linspace(0, thetamax, nt);
d = 10e-3;

nf = 3000;
v = 5850;
vs = v.*linspace(0.45, 0.75, 5);
fres = v/d/2;
f =  linspace(0, 2.5*fres, nf);

for vShear = vs
    % Calculate the reflection coefficient for all f and theta
    R = getReflectionCoeff(f, theta, v, vShear, d);
    %% Solid layer immersed in a fluid
    plotdisp(f/fres, theta, R, v/vShear)
end
end

function R = getReflectionCoeff(f, theta, v, vShear, d)
%% Material parameters
fluid1 = struct('v', 1500, 'density', 1000);
fluid3 = fluid1;
layer = struct('v', v, 'density', 7850, 'vShear', vShear);

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, d);
% Calculate V
[R, N1, N2, M1, M2] = analyticRT(f, theta, model); %#ok<ASGLU>
% Save calculated parameters
outfn = sprintf('data_%2.2f.mat', v/vShear);
save(outfn, 'f', 'theta', 'vShear', 'R', 'N1', 'N2', 'M1', 'M2'); 
end

function plotdisp(x, y, V, clcp)
%% Plot
fig = figure;
imagesc(x, y, abs(V)')
axis xy
xlabel('f/fres')
ylabel('\theta')
title(sprintf('Solid layer immersed in fluid. cL/cS = %2.2f', clcp))
colormap(gray)
outfn = fullfile('/Users/grundew/Dropbox/phd_hive/work/DispersionCurves',...
    sprintf('curves_%2.2f.fig', clcp));
saveas(fig, outfn,'fig');
end