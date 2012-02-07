function M = inputMatrix(rho, w, k_z_L)
% Function for calculating the input matrix for a fluid solid interface.
% 
% M = inputMatrix(rho, w, k_z_L)
%
% J. AcoustS. oc.A m.8 9 (4), Pt. 1, April 1991
% A new efficient algorithm to compute the exact reflection
% and transmission factors for plane waves in layered absorbing
% media (liquids and solids)
% Pierre Cervenka and Pascal Challande
% Equation 31.
a = -rho*w^2;
M = [k_z_L, -k_z_L;...
    a, a];
end