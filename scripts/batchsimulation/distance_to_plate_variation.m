%% Batch simulation of varying impdance from Z = 525 to 1.5e6
simdir = 'distanceTx_variation_aperlambda_10';
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

% Set aRx and aTx such that the fresnel distance is 1
% p.aRx = sqrt(p.cf/fres);

% Set aRx and aTx such that a/lambda = 5

p.aRx = 10*lambda;
p.aTx = p.aRx;

n = 500;
zoverfresnel = linspace(0.01, 10, n);
dist = p.aRx^2/lambda*zoverfresnel;

p.filenamevars = {'distanceTx', 'distanceRx'};

%% Do the simulations
for i = 1:n
    p.distanceTx = dist(i);
    p.distanceRx = dist(i);
    startAsmSimulation(p);
end
