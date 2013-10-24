% Script for varying cp
simdir = 'varcf-constantpoisson';
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
cpmin = 500;
cpmax = 8000;
nu = 0.2752;
cp = linspace(cpmin, cpmax, nn);
cs = cp./(sqrt(2*(1 - nu)./(1 - 2*nu)));

% Pack it up
p.filenamevars = {'cs', 'cp'};

%% Do the simulations
for i = 1:nn
    p.cp = cp(i);
    p.cs = cs(i);
    startAsmSimulation(p);
end