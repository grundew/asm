% Script for varying cf
simdir = 'varx0_water';
mkdir(simdir);
cd(simdir);

%% Get default parameters
p = parseAsmInput();
fres = 0.5*p.cp/p.thickness;
p.cf = 1500;
p.rho_fluid = 1000;
p.savemat = true;

% Define freuencies
nfft = 2^10;
fs = 2e6;
p.f = linspace(0, 1e6, nfft);

% Vary x0
x0min = 0;
x0max = 0.2;
nx0 = 20;
x0 = linspace(x0min, x0max, nx0);
nn = length(x0);

% Pack it up
p.filenamevars = {'displaceRx'};

%% Do the simulations
for i = 1:nn
    p.displaceRx = x0(i);
    startAsmSimulation(p);
end
cd('..')