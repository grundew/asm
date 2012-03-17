function [p, debug] = propagateWave(x, z, w, kx, kz, rho, V)
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

p = zeros(size(z));
q = sin(atan(kx./kz));
for i = 1:numel(z)
    % Scale Phi. Convert from velocity to pressure spectrum.
    Phi = 1i*w*rho./(1i*kz).*V;
    %Phi = k/2/pi*Phi;
    E = exp(1i*(kx.*x(i) + kz.*z(i)));
    p(i) = trapz(q, Phi.*E);
end
end