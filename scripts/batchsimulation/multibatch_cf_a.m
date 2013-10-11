% script for varying a and cf, keeping the impedance constant
simdir = 'multibatch_cf_a';
mkdir(simdir);
cd(simdir);

%% Get default parameters
p = parseAsmInput();
fres = 0.5*p.cp/p.thickness;
z = p.cf*p.rho_fluid;

nf = 2000;
f = linspace(0.7, 1.5*fres, nf);

% Define the a
na = 300;
amax = 1e-1;
amin = 1e-4;
a = logspace(log10(amin), log10(amax), na);

% Define cf and rho for constant impedance of fluid
nf = 750;
cfmin = 100;
cfmax = 3000;
cf = linspace(cfmin, cfmax, nf);
rho_fluid = z./cf;

% Pack it up
prefix = 'multibatch_cf_a';
primvarnames = {'aTx', 'aRx'};
primvar = {a, a};
secndvarnames = {'cf', 'rho_fluid'};
secndvar = {cf, rho_fluid};

%% Do the simulations
for i = 1:n
    startMultipleBatchSimulations(prefix, primvarnames, primvar,...
        secndvarnames, secndvar, p)
end
