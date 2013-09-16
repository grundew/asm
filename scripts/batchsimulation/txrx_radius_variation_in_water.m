%% Batch simulation of varying impdance from Z = 525 to 1.5e6
simdir = 'txrx_radius_variation_in_water';
mkdir(simdir);
cd(simdir);
n = 100;
radius = linspace(0.1e-3, 20e-2, n);
c = 1500;
rho = 1000;

for i = 1:n
    fnvar = {'aTx', 'aRx'};
    startAsmSimulation('cf', c,...
	    'rho_fluid', rho,...
	    'aRx', radius(i),...
	    'aTx', radius(i),...
	    'filenamevars', fnvar);
end
