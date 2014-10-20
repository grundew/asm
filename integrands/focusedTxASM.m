function X = focusedTxASM(z, k_rho, a, F, k, c, nr)
r = linspace(0, a, nr);

[KRHO, R] = meshgrid(k_rho, r);
KZ = sqrt(k^2 - KRHO.^2);
Y = 1./KZ.*R.*besselj(0, KRHO.*R).*...
        exp(1i*KZ.*(z - F + sqrt(F.^2 - R.^2)));
X = trapz(r, Y, 1);

X(k_rho==0) = 0;
X = 2*pi*c*k*X;
end