function I = integrandFluidSolidFluid_pointrx(theta_z, f, aTx,...
    c, rho, d, x0, model, reflection)
% Angular frequency and total length of wave vector
% Remember to multiply with 4*pi*k/rho/c/a
w = 2*pi*f;
k = w./c;
p = cos(theta_z);
q = sin(theta_z);
k_z = k*p;
k_rho = k*q;

%% Transmitter spatial spectrum
Tx = besselj(1, k_z*aTx);

%% Plate response, angular
% Reflection/Transmission coefficient
if reflection
    % TODO: Add loss
    Plate = analyticRTFast(w/2/pi, theta_z, model);
else
    % TODO: Add loss
    Plate = transmissionCoefficientAnalytical(f, q, model);
end

%% Phase shift from transmitter to plate and from plate to receiver
Phase_z = exp(1i*k_z*d);
F = besselj(0, k_rho*x0);

%% Integrand
I = F.*Tx.*p.^2.*Plate.*Phase_z;

if any(isnan(I))
    fprintf('NaN value detected at frequency %f and angle %f\n',...
    f, theta_z(isnan(I)));
end

end