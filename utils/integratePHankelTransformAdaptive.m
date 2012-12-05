function p = integratePHankelTransformAdaptive(ntheta, f, r, z,...
    model, a, qmax, alphaD, alphaT, debug)

if nargin < 8
    debug = false;
end

% Init p
nr = length(r);
p = zeros(nr, 1);

% Calculate integrand
q = linspace(0, qmax, ntheta+1);
dq = q(2) - q(1);
v_fluid = model.fluid(1).v;
rho_fluid = model.fluid(1).density;

% Vertical and horizontal wave number
kz = 2*pi*f*sqrt(1-q.^2)/v_fluid;
kr = 2*pi*f*q/v_fluid;

Phi = planePistonPressureAngularSpectrum(z, kr, kz, a, v_fluid, rho_fluid);
T = transmissionCoefficientAnalytical(f, q, model);
Ht = Phi.*T.*exp(1i*kz*z);

% Find the peaks
[qsid, qeid] = getPeaksStartStopIndex(abs(T));

% Sample more densily around those points
ndense = 0.5;
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
        for j = 1:nr
            Integrand = kr.*Ht.*besselj(0, kr.*r(j));
            p(j) = p(j) + integrate(Integrand(qsparseid), dq);
        end
    end
    
    % Dense indices for the peaks
    idstoppk = min(qeid(i+1), ntheta);
    qdense = linspace(q(idstop), q(idstoppk), ntheta*ndense);
    dqdense = qdense(2) - qdense(1);
    % Add the integration of the dense points
    kzdense = 2*pi*f*sqrt(1-qdense.^2)/v_fluid;
    krdense = 2*pi*f*qdense/v_fluid;

    Phidense = planePistonPressureAngularSpectrum(z, krdense, kzdense, a, v_fluid, rho_fluid);
    Tdense = transmissionCoefficientAnalytical(f, qdense, model);
    Htdense = Phidense.*Tdense.*exp(1i*kzdense*z);
    
    for j = 1:nr
        Idense = krdense.*Htdense.*besselj(0, krdense.*r(j));
        p(j) = p(j) + integrate(Idense, dqdense);
    end
    
end

end

function p = integrate(I, h)
p = h*(0.5*(I(1) + I(end)) + sum(I(2:end-1)));
end