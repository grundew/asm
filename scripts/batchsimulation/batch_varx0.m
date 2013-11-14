% Script for varying cf
simdir = 'varx0_water';
mkdir(simdir);
cd(simdir);

%% Get default parameters
p = parseAsmInput();
fres = 0.5*p.cp/p.thickness;
p.cf = 1500;
p.savemat = true;

% Define freuencies
nf = 4000;
p.f = linspace(0, 2e6, nf);

% Vary x0
nn = 50;
x0min = 0;
x0max = 20e-2;
x0 = linspace(x0min, x0max, nn);

% Pack it up
p.filenamevars = {'x0'};

%% Do the simulations
for i = 1:nn
    p.x0 = x0(i);
    startAsmSimulation(p);
end
cd('..')