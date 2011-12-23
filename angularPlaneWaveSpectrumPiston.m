function V = angularPlaneWaveSpectrumPiston(a, kx)
% P =  angularPlanewaveSpectrumPiston(a, verbose)
% 
% a - Transducer radius
% kx - Horizontal wave numbers
% ky - Wave numbers in the y-direction
% V - Angular plane wave velocity spectrum
v0 = 1;
x = kx*a;
J1 = besselj(1, x);
V = 2*pi*a^2*v0*J1./x;

% Insert correct limit when kx = 0
idzero = isnan(V);
if nnz(idzero) > 0
    V(idzero) = pi*a^2*v0;
end

end