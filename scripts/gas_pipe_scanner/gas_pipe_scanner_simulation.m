%% Simulation of gas pip scanner setup with an angle


fs = 2.5e6;
%% Transducer stuff
aTx = 6e-3;
dist = 0.08;

%% Material parameters
rho_fluid = 120;
v_fluid = 430;
thickness = 25e-3;

%% Angle
nalpha = 1;
alpha = 0;

for i = 1:nalpha
    fnvars = {'alpha'};
    startAsmSimulation_txangle('aTx', aTx, 'aRx', aTx,...
        'distanceTx', dist, 'distanceRx', dist, ...
        'fs', fs, 'filenamevars', fnvars);
end
