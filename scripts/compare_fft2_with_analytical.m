a = 6e-3;
D = 15e-3;
N = 2^13;

% Discretization
delta = D/N;
fmax = 1/2/delta; % Maximum spatial frequency with this discretisation
n = -N/2:N/2-1;
x = delta*n;
y = x;

% Circular transducer
[X, Y] = meshgrid(x, x);
un = zeros(size(X));
un(X.^2 + Y.^2 <= a^2) = 1;

% Compute wave numbers
fx = n/N/delta; % Sampling frequency in space
fy = n/N/delta; % Sampling frequency in space
frho = sqrt(fx.^2 + fy.^2); % Circular transducer
kx = 2*pi*fx;
ky = 2*pi*fy;
krho = sqrt(kx.^2 + ky.^2); % Cylindrical coordinates

% Angular spectrum
Un = fftshift(delta^2*fft2(un, N, N));

%% Analytical
[KX, KY] = meshgrid(kx, ky);
KRHO = sqrt(KX.^2 + KY.^2);
W = 2.*besselj(1, a*KRHO)/a./KRHO;
% lim x->0 J_1(x)/x -> 0.5 
W(KRHO==0) = 1;
An = pi*a^2*W;

%% Plot velocity field
figure
imagesc(x, y, un)

%% Plot image
figure
ax(1) = subplot(211);
imagesc(kx, ky, db(abs(Un)))
ax(2) = subplot(212);
imagesc(kx, ky, db(abs(An)))
linkaxes(ax, 'xy');

%% Plot line for kx = 0
figure
idkxz = (KX == 0);
plot(KRHO(idkxz), abs(Un(idkxz)), '.', KRHO(idkxz), abs(An(idkxz)), 'o')
legend('FFT2', 'Analytical')