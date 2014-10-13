function tests = focusedTxTest()
% TODO: Make a proper unittest
% TODO: Test width and amplitude
% tests = {@testFocused};
tests = functiontests(localfunctions);
end

function testFocused(testCase)

%% Parameters
a = 10e-3;
r0 = 60e-3;


%% Along acoustical axis
nz = 1000;
z = linspace(-r0, 1, nz);
nr = 4;
r_ = [0, 10, 40, 60, 82]*1e-3;
r = r_ - r0;
P = pressure(r, z, a, r0);

figure
plot(z+r0, abs(P))
set(gca, 'YScale', 'log');
title('Acoustical axis');
ylabel('Pressure')
xlabel('z (m)')
legend(arrayfun(@num2str, r + r0, 'uni', 0))

%% Along radius in the focal plane
z = [10, 25, 40, 60, 82]*1e-3 - r0;

nr = 1000;
r_ = linspace(0, 0.3, nr);
P_ = pressure(r_, z, a, r0);

r = [-fliplr(r_(2:end)), r_];
P = [flipud(P_(2:end, :)); P_];

figure
plot(r, abs(P))
set(gca, 'YScale', 'log');
title('Focal plane');
ylabel('Pressure')
xlabel('r (m)')
legend(arrayfun(@num2str, (z + r0)*1e3, 'uni', 0))


end

function P = pressure(r, z, a, r0)
nr_ = length(r);
nz_ = length(z);

f = 1.7e6;
w = 2*pi*f;
c = 1500;
k = w/c;
nz = 5000;

P = zeros(nr_, nz_);
for i = 1:nz_
    for j = 1:nr_
        func = @(theta_z) angularSpectrum(...
            theta_z, k, a, r0, r(j), z(i));
        P(j, i) = quadgk(func, 0, pi/2);
    end
end

end

function y = angularSpectrum(theta_z, k, a, r0, r, z)
k_rho = k*sin(theta_z);
k_z = k*cos(theta_z);
id = k_rho <= k*a/r0;
w = zeros(size(k_rho));
w(id) = 1;
y = 2*pi*a^2*w.*k.*sin(theta_z)*k.*cos(theta_z).*...
    exp(1i*k_z*z).*besselj(0, k_rho*r);
end