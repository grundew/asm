function I = integrandFluidSolidFluid_2Dangle_annularRx(theta_z, f, aRx, aTx,...
    c, rho, d, model, alpha, refl)
% Angular frequency and total length of wave vector
w = 2*pi*f;
k = w./c;

%% Transmitter spatial spectrum
xx = k*sin(theta_z)*aTx;
W = besselj(1, xx)./xx;
W(xx==0) = 0.5;
PhiTx = 2*W;

%% Receiver spatial spectrum
% See confluence page on Angular Spectrum Method for details on the
% relations between the angles
theta_rx = theta_z - 2*alpha;
q_rx = sin(theta_rx);
kr_rx = k*q_rx;
xx = kr_rx*aRx;
Wouter = besselj(1, xx)./xx;
Wouter(xx==0) = 0.5;
PhiRx = 2*(Wouter - W);

%% Plate response, angular
% See confluence page on Angular Spectrum Method for details on the
% relations between the angles
theta_plate = alpha - theta_z;

%% Reflection coefficient
% Todo: add loss
Plate = analyticRTFast(w/2/pi, theta_plate, model);

%% Phase shift from transmitter to plate and from plate to receiver
z = d + d*cos(2*alpha);
x = -d*sin(2*alpha);
Phase = exp(1i*k*(cos(theta_z)*z + sin(theta_z)*x));

I = rho*c*k*cos(theta_z).*cos(theta_z + 2*alpha).*Plate.*PhiRx.*PhiTx.*Phase;

if any(isnan(I))
    fprintf('NaN value detected at frequency %f and angle %f\n', f, theta_z(isnan(I)));
end

end