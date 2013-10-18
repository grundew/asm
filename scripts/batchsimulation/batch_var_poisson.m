% script for varying poisson ratio
simdir = 'varpoisson';
mkdir(simdir);
cd(simdir);

%% Get default parameters
p = parseAsmInput();
fres = 0.5*p.cp/p.thickness;

% Define freuencies
nf = 2000;
p.f = linspace(0.5*fres, 1.5*fres, nf);

% Vary Poisson's ratio
nn = 400;
numin = 0.1;
numax = 0.49;
nu = linspace(numin, numax, nn);
cp = p.cp;
cs = cp./(sqrt(2*(1 - nu)./(1 - 2*nu)));

% Pack it up
fnvar = {'cs', 'cp'};

%% Do the simulations
for i = 1:nn
    startAsmSimulation('cs', cs(i), 'filenamevars', fnvar);
end