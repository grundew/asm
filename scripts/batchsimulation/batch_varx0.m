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
x0max = 0.2;
dx0 = 0.01;
x0 = x0min:dx0:x0max;

% Pack it up
p.filenamevars = {'displaceRx'};

%% Do the simulations
for i = 1:nn
    p.displaceRx = x0(i);
    startAsmSimulation_withDisplacement(p);
end
cd('..')