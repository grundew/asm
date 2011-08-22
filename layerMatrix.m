function A = layerMatrix(rho, w, k_z_S, k_z_L, K, d)
% Function for calculating the solid layer matrix.
% 
% A = layerMatrix(rho, w, k_z)
%
% J. AcoustS. oc.A m.8 9 (4), Pt. 1, April 1991
% A new efficient algorithm to compute the exact reflection
% and transmission factors for plane waves in layered absorbing
% media (liquids and solids)
% Pierre Cervenka and Pascal Challande
% Equation 14.

% a(d) = MxP(d)xM-1
k_S = sqrt(rho*w^2/mu);
k_L = sqrt(rho*w^2/(lambda+2*mu));

S = K/k_S;
C2 = 1 - 2*S^2;

% Elements of P(d)
S_S = 1j*sin(k_z_S*d);
C_S = cos(k_z_S*d);
S_L = 1j*sin(k_z_L*d);
C_L = cos(k_z_S*d);

% Compressed notation pp. 1588
if k_z_S==0
    m_S = 0;
    d_S = 1j*d;
else
    m_S = k_z_S*S_S;
    d_S = S_S/k_z_S;
end

if k_z_L==0
    m_L = 0;
    d_L = 1j*d;
else
    m_L = k_z_L*S_L;
    d_L = S_L/k_z_L;
end

% Calculate matrix elements used twice or more
A13 = K*(C_S - C_L)/rho/w^2;
A22 = C2*C_L + 2*S^2*C_S;
A12 = K*C2*d_L - 2*S*m_L/k_S;
A11 = C2*C_S + S^2*C_L;
A21 = -(K*C_2*d_S - 2*S*m_L/k_S);
A31 = -rho*w^2*(C2*d_S + 4*S^2*m_L/(k_S)^2);

% Put together the whole shebang
A = [A11, A12, A13, -(K^2*d_L + m_S)/rho/w^2;...
     A21, A22, -(K^2*d_S + m_L)/rho/w^2, A13;...
     2*rho*w^2*S*C2*(C_S-C_L)/k_S, -rho*w^2*(C2^2*d_L + 4*S^2*m_S/(k_S)^2), A22, A12;...
     A31, A31, A21, A11];
end