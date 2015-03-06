function I = integrandFluidSolidFluid_planepiston(theta_z, f, k, d1, d3,...
    aTx, aRx, c_F, al_dB, refl, rho_F, rho_S, c_Lr, c_Sr, thick, x0)
% I = integrandFluidSolidFluid_planepiston(theta_z, f, k, d1, d3,...
%       aTx, aRx, c_F, al_dB, refl, rho_F, rho_S, c_Lr, c_Sr, thick, x0)
%
%
% Orofino
%
%
% Input:
% theta_z - Propagation angle (rad)
% f       - Frequency (Hz)
% k       - Wavenumber (1/m)
% d1      - Distance from transmitter to plate (m)
% d3      - Distance from receiver to plate (m)
% aTx     - Radius of transmitter (m)
% aRx     - Radius of receiver (m)
% c_F     - Speed of sound in fluid (m/s)
% al_dB   - Damping in the solid plate ($\alpha\lambda$ dB)
% refl    - Boolean. If true the reflection coefficient is used, otherwise
%           the transmission coefficient is used.
% rho_F   - Density of fluid (kg/m^3)
% rho_S   - Density of solid (kg/m^3)
% c_Lr    - Longitudenal speed of sound (m/s)
% c_Sr    - Shear speed of sound (m/s)
% thick   - Plate thickness (m)
% x0      - Displacement of receiver from acoustical axis (m)
%
%
% Output:
% I       - Integrand same size as theta_z


%% Compute wave numbers
q = sin(theta_z);
p = sqrt(1-q.^2);

kr = k*q;
kz = k*p;


%% Transmitter spatial sensitivity
% No angle adjustment
xx = kr*aTx;
W = besselj(1, xx)./xx;
W(xx==0) = 0.5;
PhiTx = 2*pi*aTx^2*W;


%% Receiver spatial sensitivity
xx = kr*aRx;
Wrx = besselj(1, xx)./xx;
Wrx(xx==0) = 0.5;
PhiRx = 2*pi*aRx^2*Wrx;


%% Plate response, angular
% Multiply with wave length and convert from dB to linear
% Loss parameter
if al_dB ~= 0
    % log(10)/10 = 0.2303
    alphaL = al_dB*0.2303*f/c_F;
else
    alphaL = 0;
end


%% Reflection/Transmission coefficient
if refl
    Plate = reflectionCoefficientAnalytical(f, q,...
        thick, rho_F, rho_S, c_Lr, c_Sr, c_F, alphaL);
else
    [~, Plate] = reflectionCoefficientAnalytical(f, q,...
        thick, rho_F, rho_S, c_Lr, c_Sr, c_F, alphaL);
end


%% Displacement factor
if x0 > 0
    dispRx = besselj(0, x0*kr);
else
    dispRx = 1;
end


%% Phase shift from transmitter to plate and from plate to receiver
Phase = exp(1i*kz*(d1 + d3));


%% Assemble integrand
I = Plate.*k.*q.*dispRx.*PhiRx.*PhiTx.*Phase.*k.*p.^2;

end