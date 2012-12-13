function I = hankelintegrand(q, f, r, z, model, a, alphaLambda)
v_fluid = model.fluid(1).v;
rho_fluid = model.fluid(1).density;

% Vertical and horizontal wave number
kz = 2*pi*f*sqrt(1-q.^2)/v_fluid;
kr = 2*pi*f*q/v_fluid;

Phi = planePistonPressureAngularSpectrum(z, kr, kz, a, v_fluid, rho_fluid);
alphaL = -alphaLambda*f/model.solid.v;
T = transmissionCoefficientAnalytical(f, q, model, alphaL);
Ht = Phi.*T.*exp(1i*kz*z);
I = kr.*Ht.*besselj(0, kr.*r);
end