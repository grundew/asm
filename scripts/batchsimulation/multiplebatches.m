%% Batch simulation of varying speed of sound (constant impedance) and aTx/aRx

% Vary speed of sound and keep impedance constant
nn = 20;
cf = linspace(300, 3000, nn);
z = 1000*1500;
density = z./cf;

% Tx and Rx radius
n = 100;
a = linspace(0.1e-3, 20e-2, n);

for q = 1:nn
    simdir = sprintf('s1_var_cf_%d', cf(q));
    mkdir(simdir);
    cd(simdir);
    for i = 1:n
        fnvar = {'aTx'};
        startAsmSimulation('aTx', a, 'aRx', a,...
            'rho_fluid', density(q),...
            'cf', cf(q), 'filenamevars', fnvar);
    end
end
