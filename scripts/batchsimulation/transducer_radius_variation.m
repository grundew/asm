%% Batch simulation of varying impdance from Z = 525 to 1.5e6
simdir = 'transducer_radius_variation';
mkdir(simdir);
cd(simdir);
n = 100;
radius = linspace(0.1e-3, 20e-2, n);

for i = 1:n
    fnvar = {'aTx', 'aRx'};
    startAsmSimulation('aTx', radius(i), 'filenamevars', fnvar);
end
