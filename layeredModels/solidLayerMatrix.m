function A = solidLayerMatrix(rho, w, k_z_S, k_z_L, K, k_S, d)
% A = solidLayerMatrix(rho, w, k_z_S, k_z_L, K, k_S, d)
%
%
% Function for calculating the solid layer matrix defined
% in Equation 14 in [1].
%
%
% Input:
% rho   - Density of the solid (kg/m^3)
% w     - Angular frequency (rad/s)
% k_z_S - Component of the shear wave number, z-direction (vertical, 1/m)
% k_z_L - Component of the longitudenal wave number, z-direction (1/m)
% K     - Horizontal wave number (1/m)
% k_S   - Length of shear wave number (1/m)
% d     - Thickness of layer (m)
%
%
% Output:
% A     - 4x4 Matrix
%
%
% [1] - A new efficient algorithm to compute the exact reflection
%       and transmission factors for plane waves in 
%       layered absorbing media (liquids and solids)
%       Cervenka, Pierre and Challande, Pascal,
%
%       The Journal of the Acoustical Society of America,
%       89, 1579-1589 (1991), DOI:http://dx.doi.org/10.1121/1.400993



% a(d) = MxP(d)xM-1

S = K/k_S;
C2 = 1 - 2*S^2;

% Elements of P(d), equation (9)
S_S = 1j*sin(k_z_S*d);
C_S = cos(k_z_S*d);
S_L = 1j*sin(k_z_L*d);
C_L = cos(k_z_L*d);

% Compressed notation pp. 1588
if k_z_S == 0
    m_S = 0;
    d_S = 1j*d;
else
    m_S = k_z_S*S_S;
    d_S = S_S/k_z_S;
end

if k_z_L == 0
    m_L = 0;
    d_L = 1j*d;
else
    m_L = k_z_L*S_L;
    d_L = S_L/k_z_L;
end

A = zeros(4,4);
A(1,1) = C2*C_S + 2*S^2*C_L;
A(1,2) = K*C2*d_L - 2*S*m_S/k_S;
A(1,3) = K*(C_S - C_L)/rho/w^2;
A(1,4) = -(K^2*d_L + m_S)/rho/w^2;

A(2,1) = -(K*C2*d_S - 2*S*m_L/k_S);
A(2,2) = C2*C_L + 2*S^2*C_S;
A(2,3) = -(K^2*d_S + m_L)/rho/w^2;
A(2,4) = A(1,3);

A(3,1) = 2*rho*w^2*S*C2*(C_S-C_L)/k_S;
A(3,2) = -rho*w^2*(C2^2*d_L + 4*S^2*m_S/k_S^2);
A(3,3) = A(2,2);
A(3,4) = A(1,2);

A(4,1) = -rho*w^2*(C2^2*d_S + 4*S^2*m_L/k_S^2);
A(4,2) = A(3,1);
A(4,3) = A(2,1);
A(4,4) = A(1,1);
end