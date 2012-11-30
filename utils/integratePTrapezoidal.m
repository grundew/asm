function p = integratePTrapezoidal(ntheta, f, x, z, model, a, qmax, debug)

if nargin < 8
    debug = false;
end

q = linspace(0, qmax, ntheta+1);
dq = q(2) - q(1);

v_fluid = model.fluid(1).v;
rho_fluid = model.fluid(1).density;

% Vertical and horizontal wave number
kz = 2*pi*f*sqrt(1-q.^2)/v_fluid;
kx = 2*pi*f*q/v_fluid;

% Compute the factors in the integrand
Phi = planePistonPressureAngularSpectrum(z, kx, kz, a, v_fluid, rho_fluid);
[~, T] = analyticRTFast(f, asin(q), model);
Ht = Phi.*T.*exp(1i*kz*z).*exp(1i*kx.*x);

% Trapezoidal integration
p = 2*dq*(0.5*(Ht(1) + Ht(end)) + sum(Ht(2:end-1)));

if debug
    figure
    subplot(211)
    plot(q, db(abs(Ht)), '.-')
    xlabel('q')
    ylabel('Integrand')
    tstr = sprintf('Frequency %f. N_{theta} = %d', f, ntheta);
    title(tstr)
    subplot(212)
    plot(q, unwrap(angle(Ht)), '.-')
end

end