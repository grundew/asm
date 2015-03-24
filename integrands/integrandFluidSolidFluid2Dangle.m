function I = integrandFluidSolidFluid2Dangle(theta_z, f,...
    aRx, aTx, z1, c_F, rho_F,...
    reflection, perfectReflection, al_dB, rho_S, c_Lr, c_Sr, thick, alpha)
% integrandFluidSolidFluid2Dangle - Computes the frequency response of a
% transmitter - plate - receiver, with a misalignment angle between the
% plate and the transmitter/receiver, $\alpha$.
%
%
% Note that this function should be integrated from $-\pi/2$ to $\pi/2$.
%
%
% Orofino et. al. - http://dx.doi.org/10.1121/1.405408
%
%
% Input:
% theta_z           - Propagation angle (rad)
% f                 - Frequency (Hz)
% aRx               - Radius of receiver (m). If set to zero a point
%                     source is assumed
% aTx               - Radius of transmitter (m). If set to zero a point
%                     source is assumed
% zTx               - Distance from transmitter to plate (m)
% c_F               - Speed of sound in fluid (m/s)
% rho_F             - Density of fluid (kg/m^3)
% reflection        - Boolean. If true the reflection coefficient is used,
%                     otherwise the transmission coefficient is used.
% perfectReflection - Boolean. If true, the reflection/transmission
%                     coefficient is set to 1.
% al_dB             - Damping in the solid plate ($\alpha\lambda$ dB)
% rho_S             - Density of solid (kg/m^3)
% c_Lr              - Longitudenal speed of sound (m/s)
% c_Sr              - Shear speed of sound (m/s)
% thick             - Plate thickness (m)
% alpha             - Misalignment angle between plate and transducers (rad)
%
%
% Output:
% I - Integrand, same size as theta_z


%% Angular frequency and total length of wave vector
w = 2*pi*f;
k = w./c_F;


%% Transmitter spatial spectrum
xx = k*sin(theta_z)*aTx;
W = besselj(1, xx)./xx;
W(xx==0) = 0.5;
PhiTx = 2*W;


%% Receiver spatial spectrum
% See confluence page on Angular Spectrum Method for details on the
% relations between the angles
theta_rx = theta_z + 2*alpha;
q_rx = sin(theta_rx);
kr_rx = k*q_rx;
xx = kr_rx*aRx;
W = besselj(1, xx)./xx;
W(xx==0) = 0.5;
PhiRx = 2*W;


%% Plate response, angular
% See confluence page on Angular Spectrum Method for details on the
% relations between the angles
theta_plate = alpha - theta_z;


%% Loss in plate
% Multiply with wave length and convert from dB to linear
% Loss parameter
if al_dB ~= 0
    % log(10)/10 = 0.2303
    alphaL = al_dB*0.2303*f/c_F;
else
    alphaL = 0;
end


%% Reflection/Transmission coefficient
if perfectReflection

    % If perfect reflection the reflection coefficient is just 1
    Plate = ones(size(theta_z));
    
else

    if reflection
        
        % Reflection coefficient is used
        Plate = reflectionCoefficientAnalytical(f, sin(theta_plate),...
            thick, rho_F, rho_S, c_Lr, c_Sr, c_F, alphaL);

    else
        
        % Transmission coefficient is used
        [~, Plate] = reflectionCoefficientAnalytical(f, sin(theta_plate),...
            thick, rho_F, rho_S, c_Lr, c_Sr, c_F, alphaL);
        
    end
    
end


%% Phase shift from transmitter to plate and from plate to receiver
z = z1 + z1*cos(2*alpha);
x = -z1*sin(2*alpha);
Phase = exp(1i*k*(cos(theta_z)*z + sin(theta_z)*x));

I = k*rho_F*c_F*cos(theta_z).*cos(theta_z + 2*alpha).*Plate.*PhiRx.*PhiTx.*Phase;

if any(isnan(I))
    fprintf('NaN value detected at frequency %f and angle %f\n', f, theta_z(isnan(I)));
end

end