%% Batch simulation of varying impdance from Z = 525 to 1.5e6
simdir = 'damping_variation';
mkdir(simdir);
cd(simdir);
n = 100;
alphalambda = linspace(0, 1, n);

for i = 1:n
    fnvar = {'alphaLambda_db'};
    startAsmSimulation('alphaLambda_db', alphaLambda(i), 'filenamevars', fnvar);
end
