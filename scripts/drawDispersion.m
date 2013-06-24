function drawDispersion()

nt = 4000;
thetamax = 10; % Deg
theta = linspace(0, thetamax, nt);
d = 25e-3;

nf = 6000;
v = 5850;
vs = 3200;
fres = v/d/2;

f = linspace(400e3, 800e3, nf);

for vShear = vs
    % Calculate the reflection coefficient for all f and theta
    R = getReflectionCoeff(f, theta/180*pi, v, vShear, d);
    %% Solid layer immersed in a fluid
    plotdisp(f/fres, theta, R, v/vShear)
end
end

function R = getReflectionCoeff(f, theta, v, vShear, d)
%% Material parameters
fluid1 = struct('v', 430, 'density', 120);
% fluid3 = struct('v', 2500, 'density', 1200);
fluid3 = fluid1;
layer = struct('v', v, 'density', 7850, 'vShear', vShear);

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, d);
% Calculate V
R = zeros(length(f), length(theta));
for i = 1:length(f)
    R(i, :) = fluidSolidFluidReflectionCoefficient(f(i), theta, model);
    % [~, R(i, :)] = analyticRTFast(f(i), theta, model);
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