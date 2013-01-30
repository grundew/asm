clear all;

%% Samplings stuff
nt = 1e3;
nf = 1001;
theta = linspace(0, 0.1, nt);
ff = linspace(100e3, 800e3, nf);

%% Transducer specs
aTx = 9e-3;
aRx = 3e-3;
d1 = 10e-2;
d3 = 10e-2;

%% Material parameters
rho = 1.5;
c = 350;
v_layer = 5850;
d = 10e-3;
fres = 0.5*v_layer/d;

%% Calculate the integratand
I = zeros(nt, nf);

for i = 1:nf
    freq = ff(i);
    q = sin(theta);
    p = sqrt(1-q.^2);
    w = 2*pi*freq;
    k = w./c;
    kx = k*q;
    kz = k*p;

    %% Transducer spacial sensitivities
    PhiRx = planePistonPressureAngularSpectrum(kx, aRx, c, rho);
    PhiTx = planePistonPressureAngularSpectrum(kx, aTx, c, rho);
    I(:, i) = PhiRx.*PhiTx;
end

%% Multiply with the transfer function
load('transferfunction.mat')
Htf = interp1(f, abs(Ht), ff);
HH = repmat(Htf, nt, 1);
I = HH.*I;

%% Plot the stuff
figure
imagesc(ff/fres, theta*180/pi, db(I/max(I(:))))
axis xy
colorbar
set(gca, 'CLim', [-80, 0])

% Plot contour
figure
contour(ff/fres, theta*180/pi, db(I/max(I(:))), [-30:10:0])
