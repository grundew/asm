%% Batch simulation of varying impdance from Z = 525 to 1.5e6
simdir = 'distance_to_plate_variation';
mkdir(simdir);
cd(simdir);
n = 100;
dist = linspace(0.1e-3, 20e-2, n);

for i = 1:n
    fnvar = {'distanceTx'};
    startAsmSimulation('distanceTx', dist(i), 'filenamevars', fnvar);
end
