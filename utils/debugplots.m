function debugplots(theta, freq, Phi, T, Ht, E)
% Transmission coefficient. Abs and phase
figure
ax(1) = subplot(211);
plot(theta, db(abs(T)/max(abs(T))))
title(sprintf('Transmission coefficient - f = %f', freq))
ylabel('dB')
ax(2) = subplot(212);
plot(theta, angle(T))
ylabel('angle')
xlabel('\theta')

% Incoming pressure field
figure
ax(3) = subplot(211);
plot(theta, db(abs(Phi)/max(abs(Phi))), '.')
title(sprintf('Incoming pressure field - f = %f', freq))
ylabel('dB')
ax(4) = subplot(212);
plot(theta, unwrap(angle(Phi)))
ylabel('Angle')
xlabel('\theta')

% Full integrand
figure
ax(5) = subplot(211);
plot(theta, db(abs(Ht.*E)/max(abs(Ht.*E))))
title(sprintf('Integrand - f = %f', freq))
ylabel('dB')
ax(6) = subplot(212);
plot(theta, angle(Ht.*E))
ylabel('Angle')
xlabel('\theta')

linkaxes(ax, 'x')