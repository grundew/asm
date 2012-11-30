function p = integratePHankelTransformAdaptive(ntheta, f, r, z, model, a, qmax, debug)

if nargin < 8
    debug = false;
end

% Init p
p = 0;

% Calculate integrand
q = linspace(0, qmax, ntheta+1);
dq = q(2) - q(1);
v_fluid = model.fluid(1).v;
rho_fluid = model.fluid(1).density;

% Vertical and horizontal wave number
kz = 2*pi*f*sqrt(1-q.^2)/v_fluid;
kr = 2*pi*f*q/v_fluid;

Phi = planePistonPressureAngularSpectrum(z, kr, kz, a, v_fluid, rho_fluid);
[~, T] = analyticRTFast(f, asin(q), model);
Ht = Phi.*T.*exp(1i*kz*z);
Integrand = kr.*Ht.*besselj(0, kr.*r);

if debug
    figure
    plot(q, abs(Integrand))
    hold all
end

% Find the peaks
[qsid, qeid] = getPeaksStartStopIndex(abs(T));

% Sample more densily around those points
ndense = 1;
ninc = min(40, floor(ntheta/4));
qeid = cat(2, 1, qeid+ninc);

for i = 1:length(qeid)-1
    % Non-peak indices
    idstart = qeid(i);
    idstop = max(qsid(i)-ninc, 1);
    qsparseid = idstart:idstop;
   
    if debug
        plot(q(qsparseid), abs(Integrand(qsparseid)), '.');
    end
    
    if ~isempty(qsparseid)
        p = p + integrate(Integrand(qsparseid), dq);
    end
    
    % Dense indices for the peaks
    idstoppk = min(qeid(i+1), ntheta);
    qdense = linspace(q(idstop), q(idstoppk), ntheta*ndense);
    dqdense = qdense(2) - qdense(1);
    % Add the integration of the dense points
    kz = 2*pi*f*sqrt(1-qdense.^2)/v_fluid;
    kr = 2*pi*f*qdense/v_fluid;

    Phi = planePistonPressureAngularSpectrum(z, kr, kz, a, v_fluid, rho_fluid);
    [~, T] = analyticRTFast(f, asin(qdense), model);
    Ht = Phi.*T.*exp(1i*kz*z);
    Idense = kr.*Ht.*besselj(0, kr.*r);
    p = p + integrate(Idense, dqdense);
    if debug
        plot(qdense, abs(Idense), 'o')
    end
end

end

function p = integrate(I, h)
p = h*(0.5*(I(1) + I(end)) + sum(I(2:end-1)));
end