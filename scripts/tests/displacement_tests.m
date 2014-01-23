%% Test the code for doing displacement of receiver
p = parseAsmInput();
p.f = (100e3:1e3:300e3);
p.cf = 1500;
p.rho_fluid = 1000;
p.savemat = true;
p.filenamevars = {'displaceRx'};

mkdir('disp_test')
cd('disp_test')
dis = 0:0.05:0.2;
for i = 1:length(dis);
    p.displaceRx = dis(i);
    startAsmSimulation(p);
end
cd('..')


%% Compare
