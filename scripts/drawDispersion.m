function drawDispersion()

nt = 2000;
thetamax = 0.2;
theta = linspace(0, thetamax, nt);
d = 10.15e-3;

nf = 6000;
v = 5850;
% vs = v.*linspace(0.45, 0.75, 5);
vs = 3218;
fres = v/d/2;
% f =  linspace(0.5*fres, 3*fres, nf);
f = linspace(50e3, 1e6, nf);

for vShear = vs
    % Calculate the reflection coefficient for all f and theta
    R = getReflectionCoeff(f, theta, v, vShear, d);
    %% Solid layer immersed in a fluid
    plotdisp(f/fres, theta, R, v/vShear)
end
end

function R = getReflectionCoeff(f, theta, v, vShear, d)
%% Material parameters
fluid1 = struct('v', 350, 'density', 1000);
fluid3 = fluid1;
layer = struct('v', v, 'density', 7850, 'vShear', vShear);

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, d);
% Calculate V
R = zeros(length(f), length(theta));
for i = 1:length(f)
    [~, R(i, :)] = analyticRTFast(f(i), theta, model);
end
% Save calculated parameters
% outfn = fullfile('/Users/grundew/Dropbox/phd_hive/work/DispersionCurves',...
%     sprintf('data_%2.2f.mat', v/vShear));
% save(outfn, 'f', 'theta', 'vShear',...
%     'R', 'N1', 'N2', 'M1', 'M2', 'alpha_L', 'alpha_S'); 
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