%% Script based unit test for the angular spectrum of a focused transducer
% For visual inspection only
disp('This test is for visual inspection only');

c = 1500;
f0 = 500e3;
lambda = c / f0;
omega = 2 * pi * f0;
k = 2 * pi / lambda;

F = 5 * sqrt(2) * lambda;
a = 4 * lambda;
d = sqrt(F^2 - a^2);

% Set up the coordinate grid
rmax = 2*a;
zmin = -d + F;
zmax = d + F;

dr = lambda/4;
dz = lambda/4;

ntheta = 512;
theta = linspace(0, 1, ntheta);

r = 0:dr:rmax;
z = zmin:dz:zmax;

[Z_, R_] = meshgrid(z, r);
P_ = zeros(size(Z_));

%% Compute
len = length(z);
nnn  = progress('init', 'Wait');
tic       % only if not already called
t0 = toc; % or toc if tic has already been called
tm = t0;

k = 2*pi*f0/c;
k_rho = k*sin(theta);
k_z = k*cos(theta);

nr = 200;
for ii = 1:len
    tt = ceil((toc-t0)*(len-ii)/ii);
    progress(ii/len, sprintf('%i/%i (ETA: %ds)', ii, len, tt));
    
    X = focusedTxASM(z(ii), k_rho, a, F, k, c, nr);
    % X = focusedDucer(z(ii), k_rho, a, F, k, c);
    for jj = 1:length(r)
        P_(jj, ii) = trapz(k_rho, X.*k_rho.*besselj(0, k_rho.*r(jj))); 
    end
end

%% Mirror the solution around r=0
P = [flipud(P_(r>0, :)); P_];
[Z, R] = meshgrid(z, [-fliplr(r(r>0)), r]);

%% Plot the pressure
figure
pcolor(Z*1e3, R*1e3, abs(P))
shading interp
xlabel('z (mm)')
ylabel('r (mm)')
hold on
plot(xlim(), a*ones(2, 1)*1e3, 'k--')
plot(xlim(), -a*ones(2, 1)*1e3, 'k--')
plot(F*ones(2, 1)*1e3,ylim(), 'k--')
figure
pcolor(Z*1e3, R*1e3, angle(P))
shading interp
xlabel('z (mm)')
ylabel('r (mm)')
hold on
plot(xlim(), a*ones(2, 1)*1e3, 'k--')
plot(xlim(), -a*ones(2, 1)*1e3, 'k--')
plot(F*ones(2, 1)*1e3,ylim(), 'k--')