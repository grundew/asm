%% Batch simulation of varying impdance from Z = 525 to 1.5e6
mkdir('impedance_variation');
cd('impedance_variation');
c = 342.21;
zstart = 1.5*c;
zend = 1500*1000;
n = 100;
rho = logspace(log10(zstart/c), log10(zend/c), n);

for i = 1:n
    fnvar = {'rho_fluid'};
    startAsmSimulation('rho_fluid', rho(i), 'filenamevar', fnvar);
end