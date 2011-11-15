clear
r = 0.12; % 12 cm
[X, kx, ky] = circularAperture(r, 0);
X = X/length(kx)/length(ky);
% Pick the centre line
Phi = fftshift(X(1, :));
figure
plot(kx, abs(Phi));

c = 1500;
f = 100e3;
lambda = c/f;
k = 2*pi/lambda;

q = kx/k;

% Pick out the angles less than 45 deg
id = find(abs(q) < pi/4);
q = q(id);
kx = kx(id);
Phi = Phi(id);

d = 0.01:0.01:0.1;
x = 0:0.01:0.5;

p = zeros(length(d), length(x));

for i = 1:length(d)
    
    for j = 1:length(x)
        cq = sqrt(1-q.^2);
        E = exp(1i*k*(q.*x(j) + cq*d(i)));
        p(i, :) = trapz(Phi.*E , q);
    end
end
figure
imagesc(20*log10(abs(p)))

