a = 6e-3;
fs = 2e6;
cf = 342;
D = 0.2;
f = 600e3;
lbda = cf/f;
N = 2^8;

% Solid angle for which transmitted energy includes 95% of total
phi95 = asin(2.1205*lbda/a);
theta_0 = asin(2.44*lbda/4/a);

% Discretization
% delta = 0.5*lbda*sin(phi95);
% N = 2^nextpow2((D/delta));
delta = D/N;
n = -N/2:N/2-1;
x = delta*n;
y = x;

% Circular transducer
[X, Y] = meshgrid(x, x);
un = zeros(size(X));
un(X.^2 + Y.^2 <= a^2) = 1;

% Compute wave numbers
nx = lbda*n/N/delta; % Sampling frequency in space
ny = lbda*n/N/delta; % Sampling frequency in space

% Angular spectrum
Un = delta^2*fft2(un, N, N);

%% Propagate
z = 0:0.01:1;
[NX, NY] = meshgrid(nx, ny);
NZ = fftshift(sqrt(1 - NX.^2 - NY.^2));
FZ = NZ/lbda;

idxzero = (X == 0);
uz = zeros(length(z), length(x));

for i = 1:length(z)
    uz_ = 1/delta^2*ifft2(Un.*exp(1i*2*pi*FZ*z(i)), N, N);
    % figure,pcolor(NX, NY, abs(uz_))
    
    uz(i, :) = abs(uz_(idxzero));
end

%% Plot
figure
[XZ, Z] = meshgrid(x, z);
pcolor(XZ/a, Z, db(abs(uz)))
shading interp
colormap hsv
hold all
mm = 1/tan(theta_0);
z_0 = mm*x;
plot(x/a, abs(z_0), 'k-');

%% Plot along acoustical axis
figure
plot(z, abs(uz(XZ==0)))