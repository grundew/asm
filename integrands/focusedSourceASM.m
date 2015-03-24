function X = focusedSourceASM(z, k_rho, a, r0, k, c, nr)
r = linspace(0, a, nr);

if k == 0
    X = zeros(size(k_rho));
    return;
end

[KRHO, R] = meshgrid(k_rho, r);
KZ = sqrt(k^2 - KRHO.^2);
Y = 1/1i*1./KZ.*R.*besselj(0, KRHO.*R).*...
       exp(1i*KZ.*(z - r0 + sqrt(r0.^2 - R.^2)));
X = trapz(r, Y, 1);

X(k_rho==0) = X(2);
X = 2*pi*c*k*X;
end