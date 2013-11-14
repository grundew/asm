function distance_to_plate_vary_cp_constant_nuandz(C, cp)
%% Batch simulation of varying impdance from Z = 525 to 1.5e6
simdir = sprintf('distanceTx_variation_aperlambda_%.1f', C);
mkdir(simdir);
cd(simdir);

nu = 0.2944;
cs = cp/sqrt( (2 - 2*nu) / (1 - 2*nu) );

%% Get default parameters
p = parseAsmInput();
p.thetamax = pi/2;
p.savemat = true;
p.cf = 1500;

% Get the impedance to keep constant
z = p.cp*p.rho_solid;
rho_solid = z/cp;

% Updated the speed of sound in the solid and the density
p.cp = cp;
p.cs = cs;
p.rho_solid = rho_solid;
fres = 0.5*p.cp/p.thickness;
nf = 2000;
f = linspace(0.7*fres, 1.3*fres, nf);
p.f = f;
lambda = p.cf/fres;

p.aRx = C*lambda;
p.aTx = p.aRx;

n = 500;
% zoverfresnel = logspace(-6, 0.7, n);
% dist = p.aRx^2/lambda*zoverfresnel;
dist = linspace(1e-3, 1e-2, n);
p.filenamevars = {'distanceTx', 'distanceRx'};

%% Do the simulations
for i = 1:n
    p.distanceTx = dist(i);
    p.distanceRx = dist(i);
    startAsmSimulation(p);
end

cd('..')
end