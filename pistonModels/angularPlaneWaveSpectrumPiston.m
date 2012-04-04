function V = angularPlaneWaveSpectrumPiston(a, c, q, f)
% P =  angularPlanewaveSpectrumPiston(a, verbose)
% 
% a - Transducer radius
% kx - Horizontal wave numbers
% ky - Wave numbers in the y-direction
% V - Angular plane wave velocity spectrum
k = 2*pi*f./c;
k = repmat(k(:)', length(q), 1);
q = repmat(q(:), 1, length(f));
kx = k.*q;

v0 = 1;
x = a.*kx;
J1 = besselj(1, x);
V = 2*pi*a^2*v0*J1./x;

% Insert correct limit when kx = 0
idzero = isnan(V);
if nnz(idzero) > 0
    V(idzero) = pi*a^2*v0;
end

end