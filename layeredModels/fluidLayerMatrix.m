function B = fluidLayerMatrix(rho, w, k_z_L, d)
% B = fluidLayerMatrix(rho, w, k_z_L, d)
%
%
% Function for calculating the fluid layer matrix or solid layer matrix
% when $\theta = 0$, given on page 1588 in [1].
%
% 
% Input:
% rho   - Density (kg/m^3)
% w     - Angular frequency (rad/s)
% k_z_L - Horizontal component of longitudenal wave number (1/m)
% d     - Thickness (m)
%
%
% Output:
% B     - 2x2 Matrix
%
%
% [1] - A new efficient algorithm to compute the exact reflection
%       and transmission factors for plane waves in 
%       layered absorbing media (liquids and solids)
%       Cervenka, Pierre and Challande, Pascal,
%
%       The Journal of the Acoustical Society of America,
%       89, 1579-1589 (1991), DOI:http://dx.doi.org/10.1121/1.400993

rhow2 = rho*w^2;

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