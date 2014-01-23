function [V_rx, f] = pulse_echo2D(x0)

aTx = 5e-3;
aRx = 5e-3;
d = 10e-2;
rho_f = 1000;
c_f = 1500;

nfft = 2^12;
fs = 2e6;
f = (0:nfft-1)*fs/nfft;

V_rx = zeros(size(f));

for kk = 1:numel(f)
    w = 2*pi*f(kk);
    k = w/c_f;
    func = @(k_x_) integrand(k_x_, k, d, aTx, aRx, x0);
    V_rx(kk) = quadgk(func, -inf, inf)/k;
end

V_rx = 4/rho_f/c_f*V_rx;
end

function I = integrand(k_x, k, d, aTx, aRx, x0)

k_x2 = k_x.^2;
k2 = k^2;
if k_x == 0
    S = 1;
else
    S = sin(aTx*k_x).*sin(aRx*k_x)./k_x2;
end
I = sqrt(k2 - k_x2).*S.*exp(1i*d*sqrt(k2-k_x2)).*exp(1i*k_x*x0);

end