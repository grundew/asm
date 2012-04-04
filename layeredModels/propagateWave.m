function p = propagateWave(x, z, f, q, kx, kz, rho, V)
% p = propaGatewave(x, z, w, kx, kz, rho, V)
%
% x - x-positions.
% z - z-positions, acoustical axis.
% f - Angular frequency.
% q - Sine of incoming angle.
% k - Length of wave number.
% V - Plane wave velocity spectrum.

if length(x) ~= length(z)
    error('Error:WrongInputDimensions', 'x and z must be equal length');
end
W = repmat(2*pi*f, length(q), 1);
p = zeros([size(z), length(f)]);

% Scale Phi. Convert from velocity to pressure spectrum.
Phi = 1i*rho*W./(1i*kz).*V;
%Phi = k/2/pi*Phi;
for i = 1:size(z, 1)
    for j = 1:size(z, 2)
        E = exp(1i*(kx.*x(i, j) + kz.*z(i, j)));
        p(i, j, :) = trapz(q, Phi.*E, 1);
    end
end

end