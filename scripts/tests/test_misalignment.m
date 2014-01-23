% Script for varying cf
simdir = 'varalphaplate_water';
mkdir(simdir);
cd(simdir);

%% Get default parameters
% Parameters taken from Orofino & Pedersen - IP-ASD for elastic media

dist = 7.5e-2;
nfft = 2^12;
a = 1.88e-2;

p = parseAsmInput();
p.thickness = 0.945e-2;
p.cp = 6420;
p.cs = 3040;
p.cf = 1500;
p.rho_fluid = 1000;
p.rho_solid = 2700;
p.savemat = true;
p.aRx = a;
p.aTx = a;
p.distanceRx = dist;
p.distanceTx = dist;
p.f = linspace(0, 1e6, nfft);

% Misalignment angle
alpha = (0:5)*pi/180;

% Pack it up
p.filenamevars = {'alpha_plate'};

%% Do the simulations
for i = 1:length(alpha)
    p.alpha_plate = alpha(i);
    startAsmSimulation(p);
end
cd('..')
