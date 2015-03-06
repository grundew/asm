function A = outputSolidMatrix(rho, w, k_z_S, k_z_L, K, k_S)
% A = outputSolidMatrix(rho, w, k_z_S, k_z_L, K, k_S)
%   Detailed explanation goes here

S = K/k_S;
C2 = 1 - 2*S^2;

A = [2*S/k_S, 0, -1/rho/w^2, 0;...
    0, 2*S/k_S, 0, 1/rho/w^2;...
    2*S*k_z_L/k_S, -C2, -k_z_L/rho/w^2, K/rho/w^2;...
    C2, 2*S*k_z_S/k_S, K/rho/w^2, k_z_S/rho/w^2];
end

