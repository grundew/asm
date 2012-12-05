theta_l = linspace(0, pi/2, 2^17);
rho_fluid = 1000;
c_fluid = 1500;
rho_solid = 7850;
c_l = 5850;
c_t = 3218;
[Vll, Vlt, Vtt, Vtl, Wll, Wtl, theta_fluid, theta_t]...
    = solidFluidVW(theta_l, rho_fluid, c_fluid, rho_solid, c_l, c_t);

theta = theta_l*180/pi;
figure
plot(theta, abs(Vll), theta, abs(Vlt), theta, abs(Vtt), theta, abs(Vtl))
xlabel('angle')
ylabel('abs of reflection')
title('Reflection')
legend('V_{ll}', 'V_{lt}', 'V_{tt}', 'V_{tl}')
figure
plot(theta, angle(Vll), theta, angle(Vlt), theta, angle(Vtt), theta, angle(Vtl))
xlabel('angle')
ylabel('angle of reflection')
title('Reflection')
legend('V_{ll}', 'V_{lt}', 'V_{tt}', 'V_{tl}')
figure
plot(theta, abs(Wll), theta, abs(Wtl))
xlabel('angle')
ylabel('abs of transmission')
title('Transmission')
legend('W_{ll}', 'W_{tl}')