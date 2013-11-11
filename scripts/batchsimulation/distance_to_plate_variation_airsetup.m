%% Batch simulation of distance to plate
simdir = 'distanceTx_variation_airsetup';
mkdir(simdir);
cd(simdir);

%% Get default parameters
p = parseAsmInput();
p.savemat = true;
fres = 0.5*p.cp/p.thickness;
nf = 2000;
f = linspace(0.7*fres, 1.3*fres, nf);
p.f = f;

n = 500;
dist = linspace(0.5e-2, 10e-2, n);

p.filenamevars = {'distanceTx', 'distanceRx'};

%% Do the simulations
for i = 1:n
    p.distanceTx = dist(i);
    p.distanceRx = dist(i);
    startAsmSimulation(p);
end
