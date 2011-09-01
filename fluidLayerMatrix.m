function B = fluidLayerMatrix(rho, w, k_z_L, d)
% Function for calculating the fluid layer matrix.
% 
% A = solidLayerMatrix(rho, w, k_z_S, k_z_L, K, k_S, d)
%
% J. AcoustS. oc.A m.8 9 (4), Pt. 1, April 1991
% A new efficient algorithm to compute the exact reflection
% and transmission factors for plane waves in layered absorbing
% media (liquids and solids)
% Pierre Cervenka and Pascal Challande
% Equation 29.

rhow2 = -rho*w^2;

% Elements of P(d), equation (9)
S_L = 1j*sin(k_z_L*d);
C_L = cos(k_z_L*d);

% Compressed notation pp. 1588
if k_z_L == 0
    m_L = 0;
    d_L = 1j*d;
else
    m_L = k_z_L*S_L;
    d_L = S_L/k_z_L;
end



B = [C_L, -m_L/rhow2;...
    -rhow2*d_L, C_L];
end