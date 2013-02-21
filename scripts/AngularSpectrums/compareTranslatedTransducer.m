%% Compare shifted transducer and transducer in origo
a = 6e-3;
rho = 1.5;
c = 342.21;
f = 200e3;
w = 2*pi*f;
delta_r = 10e-2;
theta = linspace(0.001, pi/2, 2^13);
kr = w/c*sin(theta);
x = kr*a;
Rx = 2*besselj(1, kr*a)./kr;
Tx = Rx*2.*besselj(0, kr*delta_r);
dr = 1/nfft/(kr(2) - kr(1));
r = (0:length(Rx)-1)*dr;


figure
subplot(211)
plot(theta*180/pi, abs(Rx), theta*180/pi, abs(Tx))
legend('Rx', 'Tx (translated)')
xlabel('\theta (degrees)')
ylabel('Angular spectrum')
xlim([0, 10])
subplot(212)
plot(r*1e3, abs(fft(Rx)), r*1e3, abs(fft(Tx)))
xlabel('r (mm)')
xlim([0 50])