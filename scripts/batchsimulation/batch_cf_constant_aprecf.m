% script for varying a and cf, keeping the impedance constant as in water
simdir = 'varcf-constant-a-per-cf';
mkdir(simdir);
cd(simdir);

%% Get default parameters
p = parseAsmInput();
fres = 0.5*p.cp/p.thickness;

% Define freuencies
nf = 1500;
p.f = linspace(0.8*fres, 1.1*fres, nf);

% Define the a
adef = 10e-3;
cdef = 1500;
aperc = adef/cdef;

% Define cf and rho for constant impedance of fluid
nc = 400;
cfmin = 100;
cfmax = 6000;
cf = linspace(cfmin, cfmax, nc);
a = aperc*cf;

% Pack it up
fnvar = {'cf', 'aTx', 'aRx'};

%% Do the simulations
for i = 1:nc
    fnvar = {'cf', 'aTx', 'aRx'};
    startAsmSimulation('cf', cf(i), 'aRx', a(i), 'aTx', a(i), 'filenamevars', fnvar);
end