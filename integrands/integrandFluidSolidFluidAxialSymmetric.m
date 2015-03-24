function I = integrandFluidSolidFluidAxialSymmetric(theta_z, f, k,...
    zTx, zRx, aTx, aRx, c_F, rho_F, displaceRx, focusRx, focusTx,...
    reflection, perfectReflection, al_dB, rho_S, c_Lr, c_Sr, thick)
% integrandFluidSolidFluidAxialSymmetric computes the integrand for an axis
% symmetric model.
%
%
%
% Orofino et. al. - http://dx.doi.org/10.1121/1.405408
%
%
% Input:
% theta_z           - Propagation angle (rad)
% f                 - Frequency (Hz)
% k                 - Wavenumber (1/m)
% zTx               - Distance from transmitter to plate (m)
% zRx               - Distance from receiver to plate (m)
% aTx               - Radius of transmitter (m). If set to zero a point
%                     source is assumed
% aRx               - Radius of receiver (m). If set to zero a point
%                     source is assumed
% c_F               - Speed of sound in fluid (m/s)
% rho_F             - Density of fluid (kg/m^3)
% displaceRx        - Displacement of receiver from acoustical axis (m)
% focusTx           - Focal raduis of transmitter. If equal to zeros,
%                     plane piston is assumed for Tx and Rx (m)
% focusRx           - Focal raduis of receiver. If equal to zeros,
%                     plane piston is assumed for Tx and Rx (m)
% reflection        - Boolean. If true the reflection coefficient is used,
%                     otherwise the transmission coefficient is used.
% perfectReflection - Boolean. If true, the reflection/transmission
%                     coefficient is set to 1.
% al_dB             - Damping in the solid plate ($\alpha\lambda$ dB)
% rho_S             - Density of solid (kg/m^3)
% c_Lr              - Longitudenal speed of sound (m/s)
% c_Sr              - Shear speed of sound (m/s)
% thick             - Plate thickness (m)
%
%
% Output:
% I - Integrand, same size as theta_z


%% Compute wave numbers
q = sin(theta_z);
p = sqrt(1-q.^2);

kr = k*q;
kz = k*p;


%% Transmitter spatial sensitivity
if aTx == 0
    
    PhiTx = ones(size(theta_z));
    
elseif aTx > 0
    
    if focusTx > 0
    
        nr = 256;
        PhiTx = focusedSourceASM(zTx, kr, aTx, focusTx, k, c_F, nr);
     
    else
    
        xx = kr*aTx;
        W = besselj(1, xx)./xx;
        W(xx==0) = 0.5;
        PhiTx = 2*pi*aTx^2*W;
            
    end
    
end


%% Receiver spatial sensitivity
if aRx == 0
    
    PhiRx = ones(size(theta_z));
    
elseif aRx > 0
    
    if focusRx > 0
    
        nr = 256;
        PhiRx = focusedSourceASM(zRx, kr, aRx, focusRx, k, c_F, nr);
     
    else
    
        xx = kr*aRx;
        W = besselj(1, xx)./xx;
        W(xx==0) = 0.5;
        PhiRx = 2*pi*aRx^2*W;
            
    end
    
end


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
        Plate = reflectionCoefficientAnalytical(f, q,...
            thick, rho_F, rho_S, c_Lr, c_Sr, c_F, alphaL);

    else
        
        % Transmission coefficient is used
        [~, Plate] = reflectionCoefficientAnalytical(f, q,...
            thick, rho_F, rho_S, c_Lr, c_Sr, c_F, alphaL);
        
    end
    
end


%% Displacement factor
if displaceRx > 0
    dispRx = besselj(0, displaceRx*kr);
else
    dispRx = 1;
end


%% Phase shift from transmitter to plate and from plate to receiver
Phase = exp(1i*kz*(zTx + zRx));


%% Assemble integrand
I = Plate.*k.*q.*dispRx.*PhiRx.*PhiTx.*Phase.*k.*p.^2;


end