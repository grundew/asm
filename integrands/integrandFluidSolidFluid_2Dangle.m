function I = integrandFluidSolidFluid_2Dangle(theta_z, f, aRx, aTx,...
    c, d1, ~, model, alpha, refl)

(theta_z, f, aRx, aTx, c, d1, d3, model, x0, alphaLambda_dB, reflection)
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
%% Reflection coefficient
% Todo: add loss
if refl
    Plate = reflectionTransmissionCoffecientAnalytical(w/2/pi, theta_plate, model);
else
    [~, Plate] = reflectionTransmissionCoffecientAnalytical(w/2/pi, theta_plate, model);
end
%% Phase shift from transmitter to plate and from plate to receiver
z = d1 + d1*cos(2*alpha);
x = -d1*sin(2*alpha);
Phase = exp(1i*k*(cos(theta_z)*z + sin(theta_z)*x));

I = k*rho*c*cos(theta_z).*cos(theta_z + 2*alpha).*Plate.*PhiRx.*PhiTx.*Phase;

if any(isnan(I))
    fprintf('NaN value detected at frequency %f and angle %f\n', f, theta_z(isnan(I)));
end

end