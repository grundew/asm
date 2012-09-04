function p = propagateReflectedWave(x, z, f, q, kx, kz, rho, V, R)
% p = propaGatewave(x, z, w, kx, kz, rho, V)
%
% x - x-positions.
% z - z-positions, acoustical axis.
% w - Angular frequency.
% q - Sine of incoming angle.
% kz - Vertical wave number
% kx - Horizontal wave number
% V - Plane wave velocity spectrum.
% rho - Density of propagation medium
% R - Transmission/reflection coefficient

% TODO: Handle q = -1 and 1 (Phi goes now to inf)
if length(x) ~= length(z)
    error('Error:WrongInputDimensions', 'x and z must be equal length');
end

% W = repmat(2*pi*f, length(q), 1);
p = zeros([size(z), length(f)]);

nq = length(q);

% Scale Phi. Convert from velocity to pressure spectrum.
% if nq == 1 || isequal(kz, zeros(size(kz)))
%     Phi = ones(size(V));
% else
%     Phi = 1i*rho*W./(1i*kz).*V;
% end
% Phi = k/2/pi*Phi;
Phi = V;

for i = 1:size(z, 1)
    for j = 1:size(z, 2)
        E = exp(1i*(kx.*x(i, j) + kz.*z(i, j)));
        if nq == 1
            p(i, j, :) = R.*Phi.*E;
        else
            p(i, j, :) = trapz(q, R.*Phi.*E, 1);
        end
    end
end

end