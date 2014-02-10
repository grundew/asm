function I = integrandFluidSolidFluid_2Dangle(theta_z, f, aRx, aTx,...
    c, rho, d, model, alpha, reflection)
% Angular frequency and total length of wave vector
w = 2*pi*f;
k = w./c;

%% Indices of angles less than alpha (the misaligment angle of the plate)
idltalpha = theta_z <= alpha;

%% Transmitter spatial spectrum
xx = k*sin(theta_z)*aTx;
W = besselj(1, xx)./xx;
W(xx==0) = 0.5;
PhiTx = 2*W;

%% Receiver spatial spectrum
% See confluence page on Angular Spectrum Method for details on the
% relations between the angles
theta_rx = theta_z - 2*alpha;
theta_rx(idltalpha) = 2*alpha - theta_z(idltalpha);
q_rx = sin(theta_rx);
kr_rx = k*q_rx;

xx = kr_rx*aRx;
W = besselj(1, xx)./xx;
W(xx==0) = 0.5;
PhiRx = 2*W;

%% Plate response, angular
% See confluence page on Angular Spectrum Method for details on the
% relations between the angles
theta_plate = theta_z - alpha;
%% Reflection/Transmission coefficient
if reflection
    % TODO: Add loss
    Plate = analyticRTFast(w/2/pi, theta_plate, model);
else
    % TODO: Add loss
    alphaL = 0;
    Plate = transmissionCoefficientAnalytical(f, sin(theta_plate), model, alphaL);
end
%% Phase shift from transmitter to plate and from plate to receiver
z = d + d*cos(2*alpha);
x = -d*sin(2*alpha);
Phase = exp(1i*k*(cos(theta_z)*z + sin(theta_z)*x));

I = k/rho/c*Plate.*PhiRx.*PhiTx.*Phase.*k.^2.*cos(theta_z);

if any(isnan(I))
    fprintf('NaN value detected at frequency %f and angle %f\n', f, theta_z(isnan(I)));
end

end