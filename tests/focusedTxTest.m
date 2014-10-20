function tests = focusedTxTest()
% TODO: Test width and amplitude
% tests = {@testFocused};
% tests = functiontests(localfunctions);

testFocused();
end

function testFocused(testCase)

%% Parameters
a = 9e-3;
r0 = 100e-3;
f = 100e3;
c = 342;


%% Along acoustical axis
nz = 1000;
z = linspace(-r0, 1, nz);
r = linspace(0, 0.1, 10);
P = pressure(r, z, a, r0, f, c);

figure
plot(z+r0, abs(P))
hold on
set(gca, 'YScale', 'log');
yl = ylim();
plot(ones(2, 1)*r0, yl, 'k')

title('Acoustical axis');
ylabel('Pressure');
xlabel('Distance from transmitter (m)');
legend(arrayfun(@num2str, r*1e3, 'uni', 0));

%% Along radius in the focal plane
z = linspace(0, 200, 5)*1e-3;

nr = 1000;
r_ = linspace(0, 0.3, nr);
P_ = pressure(r_, z, a, r0, f, c);

r = [-fliplr(r_(2:end)), r_];
P = [flipud(P_(2:end, :)); P_];

figure
plot(r, abs(P))
set(gca, 'YScale', 'log');
title('r-plane');
ylabel('Pressure')
xlabel('r (m)')
legend(arrayfun(@num2str, (z + r0)*1e3, 'uni', 0))

%% Plot time signals
z = 100e-3;
nr = 5;
r = linspace(0, 100e-3, nr);
f0 = 100e3;
fs = 1e6;
t = 0:1/fs:5/f0;
nfft = 2^12;
f_compute = (0:nfft/2+1)*fs/nfft;
x_ex = tukeywin(length(t), 0.5)'.*sin(2*pi*f0*t);
X_ex = ifft(conj(hilbert(x_ex)), nfft);

x = zeros(length(r), nfft);
for i = 1:length(r)
    P = zeros(1, nfft);
    for j = 1:length(f_compute)
        P(j) = pressure(r(i), z, a, r0, f_compute(j), c);
    end
    % Convolve with excitation pulse
    x(i, :) = real(fft(P.*X_ex, nfft));
end

figure
plot((0:size(x, 2)-1)/fs, x)

end


function P = pressure(r, z, a, r0, f, c)
nr = length(r);
nz = length(z);

w = 2*pi*f;
k = w/c;

P = zeros(nr, nz);
for i = 1:nz
    for j = 1:nr
        func = @(theta_z) angularSpectrum(...
            theta_z, k, a, r0, r(j), z(i));
        P(j, i) = quadgk(func, 0, pi/2);
    end
end

end


function y = angularSpectrum(theta_z, k, a, r0, r, z)
k_rho = k*sin(theta_z);
k_z = k*cos(theta_z);

xx = k_rho*a;
W = besselj(1, xx)./xx;
W(xx==0) = 0.5;
w = 2*pi*a^2*W;


% id = k_rho <= k*a/r0;
% w = zeros(size(k_rho));
% w(id) = 1;
y = 2*pi*a^2*w.*k.*sin(theta_z)*k.*cos(theta_z).*...
    exp(1i*k_z*z).*exp(1i*k_rho*r).*besselj(0, k_rho*r);
end