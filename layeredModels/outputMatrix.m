function M = outputMatrix(rho, w, k_z_L)
% Function for calculating the solid layer matrix.
% 
% M = outputMatrix(rho, w)
%
% J. AcoustS. oc.A m.8 9 (4), Pt. 1, April 1991
% A new efficient algorithm to compute the exact reflection
% and transmission factors for plane waves in layered absorbing
% media (liquids and solids)
% Pierre Cervenka and Pascal Challande
% Equation 33.
% V'x(N^(n+1))^-1
rhow2 = rho*w^2;
M = [0, -1/rhow2;...
    -1, -k_z_L/rhow2];
end