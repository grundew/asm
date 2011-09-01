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
ylabel('Horizontal wave number')