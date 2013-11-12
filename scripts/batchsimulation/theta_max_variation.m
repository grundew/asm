%% Batch simulation of varying impdance from Z = 525 to 1.5e6
simdir = 'thetamax_variation_aperlambda_4';
mkdir(simdir);
cd(simdir);

%% Get default parameters
p = parseAsmInput();
p.savemat = true;
p.cf = 1500;
fres = 0.5*p.cp/p.thickness;
nf = 2000;
f = linspace(0.7*fres, 1.3*fres, nf);
p.f = f;
lambda = p.cf/fres;

% Set aRx and aTx such that a/lambda = 4
p.aRx = 4*lambda;
p.aTx = p.aRx;

zoverfresnel = 0.5;
dist = p.aRx^2/lambda*zoverfresnel;
p.distanceTx = dist;
p.distanceRx = dist;

% Vary integration limit
% Critical angle, shear is 0.48, comp 0.25
n = 20;
thetamax = linspace(0.5, pi/2, n);
p.filenamevars = {'thetamax'};

%% Do the simulations
for i = 1:n
    startAsmSimulation(p);
end
