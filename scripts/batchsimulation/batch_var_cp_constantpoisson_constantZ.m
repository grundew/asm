% script for varying cp, keeping the impedance and poisson ratio constant
simdir = 'varcp_constantmu_constantz';
mkdir(simdir);
cd(simdir);

%% Get default parameters
p = parseAsmInput();
fres = 0.5*p.cp/p.thickness;

% Define freuencies
nf = 4000;
p.f = linspace(0.5*fres, 1.5*fres, nf);
z = p.cp*p.rho_solid;

% Vary cp
nn = 400;
cpmin = 500;
cpmax = 8000;
nu = 0.2752;
cp = linspace(cpmin, cpmax, nn);

% Keep poisson's ratio constant
cs = cp./(sqrt(2*(1 - nu)./(1 - 2*nu)));
% Keep impedance constant
rho = z./cp;

% Variables written to file
p.filenamevars = {'cp', 'cs', 'rho_solid'};

%% Do the simulations
for i = 1:nn
    p.cp = cp(i);
    p.cs = cs(i);
    p.rho_solid = rho(i);
    startAsmSimulation(p);
end