% script for varying a and cf, keeping the impedance constant as in water
simdir = 'multibatch-cf-a';
mkdir(simdir);
cd(simdir);

%% Get default parameters
p = parseAsmInput();
fres = 0.5*p.cp/p.thickness;
z = 1000*1500;

nf = 1500;
p.f = linspace(0.8*fres, 1.1*fres, nf);

% Define the a
na = 100;
logamax = 0;
logamin = -4;
a = logspace(logamin, logamax, na);

% Define cf and rho for constant impedance of fluid
nf = 100;
cfmin = 100;
cfmax = 6000;
cf = linspace(cfmin, cfmax, nf);
rho_fluid = z./cf;

% Pack it up
prefix = 'multibatch_cf_a';
primvarnames = {'aTx', 'aRx'};
primvar = {a, a};
secndvarnames = {'cf', 'rho_fluid'};
secndvar = {cf, rho_fluid};

%% Do the simulations
startMultipleBatchSimulations(prefix, primvarnames, primvar,...
        secndvarnames, secndvar, p)