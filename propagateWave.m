function [p, debug] = propagateWave(z, x, w, kx, kz, rho, V)
% p = propaGatewave(x, z, w, kx, kz, rho, V)
%
% x - x-positions.
% z - z-positions, acoustical axis.
% w - Angular frequency.
% q - Sine of incoming angle.
% k - Length of wave number.
% V - Plane wave velocity spectrum.

if length(x) ~= length(z)
    error('Error:WrongInputDimensions', 'x and z must be equal length');
end

p = zeros(length(z), length(x));
q = sin(atan(kx./kz));
for i = 1:length(z)
    for j = 1:length(x)
        % kx = k*q
        % kz = k*cq
        % fprintf('%d\t %d\n', z(i), x(j))
        % Scale Phi. Convert from velocity to pressure spectrum.
        Phi = 1i*w*rho./(1i*kz).*V;
        %Phi = k/2/pi*Phi;
        E = exp(1i*(kx.*x(j) + kz.*z(i)));
        p(i, j) = trapz(q, Phi.*E);
    end
end
end