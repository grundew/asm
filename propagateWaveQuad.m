function p = propagateWaveQuad(z, x, w, k, q, rho, a)
% p = propaGatewave(x, z, V)
%
% x - x-positions.
% z - z-positions, acoustical axis.
% w - Angular frequency.
% q - Sine of incoming angle.
% k - Length of wave number.
% V - Plane wave velocity spectrum.
p = zeros(length(z), length(x));

for i = 1:length(z)
    for j = 1:length(x)
        % kx = k*q
        % kz = k*cq
        % fprintf('%d\t %d\n', z(i), x(j))
        % Scale Phi. Convert from velocity to pressure spectrum.
        % Phi = w*rho./kz.*V;
        % Phi = k/2/pi*Phi;
        xx = x(j);
        zz = z(i);
        %E = exp(1i*k*(q.*x(j) + cq*z(i)));
        %p(i, j) = trapz(q, Phi.*E);
        p(i, j) = quadgk(@(x) (integrand(x)), min(q), max(q));
    end
end

    function y = integrand(sintheta)
        costheta = sqrt(1 - sintheta.^2);
        kz = k*costheta;
        Phi = angularPlaneWaveSpectrumPiston(a, k*sintheta);
        Phi = w*rho./kz.*Phi;
        Phi = k/2/pi*Phi;
        E = exp(1i*k*(sintheta.*xx + costheta.*zz));
        y = Phi.*E;
    end
end