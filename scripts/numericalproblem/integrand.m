function I = integrand(q, f, a, x, z, model)
v_fluid = model.fluid(1).v;
rho_fluid = model.fluid(1).density;
kz = 2*pi*f*sqrt(1-q.^2)/v_fluid;
kx = 2*pi*f*q/v_fluid;
Phi = planePistonPressureAngularSpectrum(z, kx, kz, a, v_fluid, rho_fluid);
[~, T] = analyticRTFast(f, asin(q), model);
I = Phi.*T.*exp(1i*kz*z).*exp(1i*kx.*x);
end