function B = fluidSolidFluidLayerMatrix(rho, w, k_z_S, k_z_L, K, k_S, d, thresh)
% Function for calculating the solid layer matrix.
%
% B = fluidSolidFluidLayer(rho, w, k_z_S, k_z_L, K, d, mu)
%
% J. AcoustS. oc.A m.8 9 (4), Pt. 1, April 1991
% A new efficient algorithm to compute the exact reflection
% and transmission factors for plane waves in layered absorbing
% media (liquids and solids)
% Pierre Cervenka and Pascal Challande
% Equation 24 pp. 1582

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

% Check if there is an evanescent longitudenal wave
% k_z_L = k_r_L + jalpha_L, exp(-alpha_L*D)<<1
alpha_L = imag(k_z_L);
alpha_S = imag(k_z_S);

if exp(-alpha_L*d) < thresh && exp(-alpha_S*d) > thresh
    % This propagates the shear wave only
    rhow2 = rho*w^2;
    
    B1 = [2*S^2 , 0;...
        rhow2*C2^2/k_z_L, 1];
    
    b = [C_S, -k_S^2*d_S/(2*rhow2);...
        -2*rhow2*m_S/k_S^2, C_S];
    
    B0 = [0.5*1/S^2, 0;...
        rhow2*C2^2/(2*S^2*k_z_L), 1];
    
    B = B0*b*B1;
    return
elseif exp(-alpha_S*d) < thresh
    B = zeros(2, 2);
    return
end

A21 = -(K*C2*d_S - 2*S*m_L/k_S);
A22 = C2*C_L + 2*S^2*C_S;
A23 = -(K^2*d_S + m_L)/(rho*w^2);

A31 = 2*rho*w^2*S*C2*(C_S - C_L)/k_S;
A32 = -rho*w^2*(C2^2*d_L + 4*S^2*m_S/k_S^2);

A41 = -rho*w^2*(C2^2*d_S + 4*S^2*m_L/k_S^2);
A42 = A31;
A43 = A21;

B(1, 1) = A22 - A21*A42/A41;
B(1, 2) = A23 - A21*A43/A41;
B(2, 1) = A32 - A31*A42/A41;
B(2, 2) = A22 - A31*A43/A41;
end