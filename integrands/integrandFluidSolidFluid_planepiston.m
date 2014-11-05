function I = integrandFluidSolidFluid_planepiston(theta_z, f, k, d1, d3, aTx, aRx, c_F,...
    al_dB, refl, rho_F, rho_S, c_Lr, c_Sr, thick, x0)
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