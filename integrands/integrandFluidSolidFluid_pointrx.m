function I = integrandFluidSolidFluid_pointrx(theta_z, f, aTx,...
    c, rho, d, x0, model, reflection)
% Angular frequency and total length of wave vector
w = 2*pi*f;
k = w./c;
q = sin(theta_z);
p = cos(theta_z);

%% Transmitter spatial spectrum
xx = k*q*aTx;
W = besselj(1, xx)./xx;
W(xx==0) = 0.5;
PhiTx = 2*W;

%% Plate response, angular
%% Reflection/Transmission coefficient
if reflection
    % TODO: Add loss
    Plate = analyticRTFast(w/2/pi, theta_z, model);
else
    % TODO: Add loss
    Plate = transmissionCoefficientAnalytical(f, q, model);
end

%% Phase shift from transmitter to plate and from plate to receiver
Phase = exp(1i*k*(p*d + q*x0));

I = Plate.*k.*q.*p/rho/c.*PhiTx.*Phase*k.*p;

if any(isnan(I))
    fprintf('NaN value detected at frequency %f and angle %f\n', f, theta_z(isnan(I)));
end

end