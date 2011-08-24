function B = fluidLayerMatrix(rho, w, k_z_L)
% Function for calculating the fluid layer matrix.
% 
% A = solidLayerMatrix(rho, w, k_z_S, k_z_L, K, k_S, d)
%
% J. AcoustS. oc.A m.8 9 (4), Pt. 1, April 1991
% A new efficient algorithm to compute the exact reflection
% and transmission factors for plane waves in layered absorbing
% media (liquids and solids)
% Pierre Cervenka and Pascal Challande
% Equation 31.

rhow2 = -rho*w^2;
B = [k_z_L, k_z_L;...
    rhow2, rhow2];
end