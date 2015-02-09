function I = integrandFluidSolidFluid_pointrx(theta_z, f, aTx,...
    c, rho, d, x0, model, alpha, reflection, al_dB)
% Angular frequency and total length of wave vector
% Remember to multiply with 4*pi*k/rho/c/a
w = 2*pi*f;
k = w./c;
p = cos(theta_z);
q = sin(theta_z);
k_z = k*p;
k_rho = k*q;


%% Loss parameter
if al_dB ~= 0
    % log(10)/10 = 0.2303
    alphaL = al_dB*0.2303*f/c_F;
else
    alphaL = 0;
end


%% Transmitter spatial spectrum
xx = k*sin(theta_z+alpha)*aTx;
W = besselj(1, xx)./xx;
W(xx==0) = 0.5;
Tx = 2*pi*aTx^2*W;


%% Plate response, angular
% Reflection/Transmission coefficient
if reflection
    % TODO: Add loss
    Plate = analyticRTFast(w/2/pi, theta_z, model);
else
    % TODO: Add loss
    Plate = transmissionCoefficientAnalytical(f, q, model, alphaL);
end


%% Phase shift from transmitter to plate and from plate to receiver
Phase_z = exp(1i*k_z*d);
F = 2*pi*besselj(0, k_rho*x0);


%% Integrand
I = k.^2*F.*Tx.*p.^2.*q.*Plate.*Phase_z/rho/c;

if any(isnan(I))
    fprintf('NaN value detected at frequency %f and angle %f\n',...
        f, theta_z(isnan(I)));
end

end