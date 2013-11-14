function distance_to_plate_variation(C)
%% Batch simulation of varying impdance from Z = 525 to 1.5e6
simdir = sprintf('distanceTx_variation_aperlambda_%.1f', C);
mkdir(simdir);
cd(simdir);

%% Get default parameters
p = parseAsmInput();
p.thetamax = pi/2;
p.savemat = true;
p.cf = 1500;
fres = 0.5*p.cp/p.thickness;
nf = 2000;
f = linspace(0.7*fres, 1.3*fres, nf);
p.f = f;
lambda = p.cf/fres;

p.aRx = C*lambda;
p.aTx = p.aRx;

n = 550;
zoverfresnel = logspace(-6, 0.7, n);
dist = p.aRx^2/lambda*zoverfresnel;
p.filenamevars = {'distanceTx', 'distanceRx'};

%% Do the simulations
for i = 1:n
    p.distanceTx = dist(i);
    p.distanceRx = dist(i);
    startAsmSimulation(p);
end

end
