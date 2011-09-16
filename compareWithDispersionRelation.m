%% Clear stuff
clear all;

%% Reflection coefficient
d = 12.4e-3;
fluid = materials.Fluid;
fluid.v = 1450;
fluid.density = 1000;
solid = materials.MaterialFactory.produce('stainless steel');
model = MultiLayerModel(fluid, solid, fluid, d);

f = 1:0.001e3:1000e3;
theta = 1e-3;
[V, W] = fluidSolidFluid(f, theta, model);
figure
subplot(211)
plot(2*pi*f, real(V));
subplot(212)
plot(2*pi*f, imag(V));

lpo = LocalPeaksPlain();
lpo.findPeaks(imag(V));
subplot(212)
hold all
plot(2*pi*f(lpo.k), lpo.x, 'x')

%% Dispersion curves
pdr = dispersion.PlateDispersionRelation(...
    solid, d);
% Horizontal wavenumber
%k_hor = linspace(0.1, 15*pi/pdr.thickness, 150);
k_hor = linspace(0.01, 15*pi/pdr.thickness, 150);
fd = linspace(1, 1200e3, 5000);
[ba, bs] = pdr.branches(fd, k_hor);

figure
hold all
for i = 1:length(ba)
    ha = plot(ba(i).frequency, ba(i).waveNumber);
end
for i = 1:length(bs)
    hs = plot(bs(i).frequency, bs(i).waveNumber, '--');
end
legend([ha, hs], 'Antisymmetric', 'Symmetric');

fs = zeros(length(bs)-1, 1);
for i = 2:length(bs)
    k_z = 2*pi*bs(i).frequency(1)/fluid.v*sin(theta);
    fs(i) = bs(i).freq(k_z);
end

fa = zeros(length(ba)-1, 1);
for i = 2:length(ba)
    k_z = 2*pi*ba(i).frequency(1)/fluid.v*sin(theta);
    fa(i) = ba(i).freq(k_z);
end
subplot(212)
hold all
plot(fa, zeros(size(fa)), 'x')
plot(fs, zeros(size(fs)), 'o')

subplot(211)
hold all
plot(fa, ones(size(fa)), 'x')
plot(fs, ones(size(fs)), 'o')


clear all;
d = 12.4e-3;
pdr = dispersion.PlateDispersionRelation(materials.MaterialFactory.produce('stainless steel'), d);
% Horizontal wavenumber
k_hor = linspace(0.1, 15*pi/pdr.thickness, 150);
f = linspace(1, 1200e3, 5000);
[ba, bs] = pdr.branches(f, k_hor);

% Total wavenumber
k = 2*pi*f/1500;
theta = linspace(0, pi/2, 100);

[R, T, KK] = waterSteelWaterC1(f, theta, d);

RR = angleToWaveNumber(f, theta, R, k_hor);

%% Plot the all the stuff
fig = figure;
imagesc(f, k_hor, abs(RR)')
axis('xy')
colormap('jet')
xlabel('Frequency (Hz)')
ylabel('Horizontal wave number')
hold all
for i = 1:length(ba)
    plot(ba(i).frequency/2/pi, ba(i).waveNumber, 'k.', 'MarkerSize', 5)
end
xlabel('Frequency (Hz)')

