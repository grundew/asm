theta = linspace(0, pi/2, 2^17);
rho_fluid = 1000;
c_fluid = 1500;
rho_solid = 7850;
c_l = 5850;
c_t = 3218;
[V, Wl, Wt] = fluidSolidVW(theta, rho_fluid, c_fluid, rho_solid, c_l, c_t);
plot(theta, abs(V), theta, abs(Wl), theta, abs(Wt))
figure
plot(theta, imag(V), theta, imag(Wl), theta, imag(Wt))