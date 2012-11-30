function p = integratePTrapezoidalStepWise(ntheta, f, x, z, model, a, qmax, nstep, debug)

if nargin < 9
    debug = false;
end

qq = linspace(0, qmax, ntheta+1);
dq = qq(2) - qq(1);

v_fluid = model.fluid(1).v;
rho_fluid = model.fluid(1).density;

if debug
    figure
    ax(1) = subplot(211);
    hold(ax(1), 'all');
    xlabel('q')
    ylabel('Integrand')
    tstr = sprintf('Frequency %f. N_{theta} = %d', f, ntheta);
    title(tstr)
    ax(2) = subplot(212);
    hold(ax(2), 'all');
    linkaxes(ax, 'x')
end


for j = 1:nstep+1
    % Add half of the sum of the end points if j == 1
    % Else add the sum of the interior points
    if j == 1
        q = [0 qmax];
    else
        q = qq(j:nstep:end-1);
    end
    % Vertical and horizontal wave number
    kz = 2*pi*f*sqrt(1-q.^2)/v_fluid;
    kx = 2*pi*f*q/v_fluid;
    
    % Compute the factors in the integrand
    Phi = planePistonPressureAngularSpectrum(z, kx, kz, a, v_fluid, rho_fluid);
    [~, T] = analyticRTFast(f, asin(q), model);
    Ht = Phi.*T.*exp(1i*kz*z).*exp(1i*kx.*x);

    % Add half the end points or the sum of the interior points
    if j == 1
        I = 0.5*sum(Ht);
    else
        I = I + sum(Ht);
    end
    
    if debug
        plot(ax(1), q, db(abs(Ht)), '.')
        plot(ax(2), q, unwrap(angle(Ht)), '.')
    end
    
end

p = 2*dq*I;
end