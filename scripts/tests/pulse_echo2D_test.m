function [V_rx, VV_rx ,f] = pulse_echo2D_test(x0)

aTx = 5e-3;
aRx = 5e-3;
d = 10e-2;
dx = 0.5e-3;
rho_f = 1000;
c_f = 1500;

f = linspace(100e3, 500e3, 2^8);
V_rx = zeros(size(f));

% for kk = 1:numel(f)
%     w = 2*pi*f(kk);
%     k = w/c_f;
% 
%     x = x0-aRx:dx:x0+aRx;
%     I = zeros(size(x));
%     for i = 1:numel(x)
%         x_ = x(i);
%         func = @(k_x_) exp(1i*k_x_*x_).*sqrt(k^2 - k_x_.^2).*sin(aTx*k_x_)./k_x_.*exp(1i*d*sqrt(k^2-k_x_.^2));
%         I(i) = 2/rho_f/c_f/k*quadgk(func, -inf, inf, 'MaxIntervalCount', 1280);
%     end
%     
%     V_rx(kk) = trapz(x, I);
%     
% end

VV_rx = zeros(size(f));

for kk = 1:numel(f)
    w = 2*pi*f(kk);
    k = w/c_f;
    func = @(k_x_) sqrt(k^2 - k_x_.^2).*sin(aTx*k_x_)./k_x_.*exp(1i*d*sqrt(k^2-k_x_.^2)).*sin(aRx*k_x_)./k_x_.*exp(1i*k_x_*x0);
    VV_rx(kk) = 4/rho_f/c_f/k*quadgk(func, -inf, inf, 'MaxIntervalCount', 1280);
end

figure
subplot(211)
plot(f*1e-3, abs(V_rx), f*1e-3, abs(VV_rx))
subplot(212)
plot(f*1e-3, unwrap(angle(V_rx)), f*1e-3, unwrap(angle(VV_rx)))
end