% script for varying poisson ratio
simdir = 'varcf';
mkdir(simdir);
cd(simdir);

%% Get default parameters
p = parseAsmInput();
fres = 0.5*p.cp/p.thickness;

% Define freuencies
nf = 4000;
p.f = linspace(0.5*fres, 1.5*fres, nf);
z = p.cp*p.rho_solid;

% Vary Poisson's ratio
nn = 1000;
cfmin = 50;
cfmax = 8000;
cf = linspace(cfmin, cfmax, nn);

% Pack it up
fnvar = {'cf'};

%% Do the simulations
for i = 1:nn
    startAsmSimulation('cf', cf(i), 'filenamevars', fnvar);
end