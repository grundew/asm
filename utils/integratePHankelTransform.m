function p = integratePHankelTransform(ntheta, f, r, z, model, a, qmax, debug)

if nargin < 8
    debug = false;
end

q = linspace(0, qmax, ntheta+1);
dq = q(2) - q(1);

v_fluid = model.fluid(1).v;
rho_fluid = model.fluid(1).density;

% Vertical and horizontal wave number
kz = 2*pi*f*sqrt(1-q.^2)/v_fluid;
kr = 2*pi*f*q/v_fluid;

% Compute the factors in the integrand
Phi = planePistonPressureAngularSpectrum(z, kr, kz, a, v_fluid, rho_fluid);
[~, T] = analyticRTFast(f, asin(q), model);
Ht = Phi.*T.*exp(1i*kz*z);
Integrand = kr.*Ht.*besselj(0, kr.*r);

p = dq*(0.5*(Integrand(1) + Integrand(end)) + sum(Integrand(2:end-1)));
if debug
    figure
    ax(1) = subplot(211);
    plot(q, abs(Integrand), '.-')
    ax(2) = subplot(212);
    plot(q, unwrap(angle(Integrand)), '.-')
    linkaxes(ax, 'x')
end
end