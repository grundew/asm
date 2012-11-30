%% Trackdown of the numerical problem

%% Set the material parameters
cF = 1500;
rho_fluid = 1000;
rho_layer = 7850;
cL = 5850;
cS = 3218;
fluid1 = struct('v', cF, 'density', rho_fluid);
fluid3 = fluid1;
layer = struct('v', cL, 'density', rho_layer, 'vShear', cS);
d = 10.15e-3;
fres = 0.5*cL/10.15e-3;
thresh = 1e-25;
ntheta = 2^12;
theta = linspace(0.001, pi/2-0.01, ntheta);
f = 2000e3;

%% Calculate the reflection and transmission coefficient
model = MultiLayerModel(fluid1, layer, fluid3, d);
[R, T] = fluidSolidFluidReflectionCoefficient(f, theta, model, thresh);

%% Calculate wave numbers and angles in the steel
w = 2*pi*f;
kF = w/cF;
kL = w/cL;
kS = w/cS;
thetaL = asin(cL/cF*sin(theta));
thetaS = asin(cS/cF*sin(theta));
kFz = kF*cos(theta);
kLz = kL*cos(thetaL);
kSz = kS*cos(thetaS);
K = kL*sin(thetaL);

% Convenience variables
S = K/kS;
C2 = 1 - 2*S.^2;
Cl = cos(kLz*d);
Cs = cos(kSz*d);
Sl = 1i*sin(kLz*d);
Ss = 1i*sin(kSz*d);
dS = Ss./kSz;
dL = Sl./kLz;
mL = kLz.*Sl;
mS = kSz.*Ss;

%% Calculate the relevant matrix elements of the A-matrix
a22 = C2.*Cl + 2*S.^2.*Cs;
a21 = -(K.*C2.*dS - 2*S.*mL./kS);
a31 = 2*rho_layer*w^2.*S.*C2.*(Cs-Cl)./kS;
a41 = -rho_layer*w^2*(C2.^2.*dS + 4*S.^2.*mL./kS.^2);
a33 = a22;
a42 = a31;
a43 = a21;
a23 = -(K.^2.*dS+mL)/rho_layer/w^2;
a32 = -rho_layer*w^2*(C2.^2.*dL + 4*S.^2.*mS./kS.^2);
b11 = a22 - a21.*a42./a41;
b22 = a33 - a31.*a43./a41;
b12 = a23 - a21.*a43./a41;
b21 = a32 - a31.*a42./a41;

%% Reflection coefficient calculated from B
alpha = (b12 + kFz./rho_fluid/w.^2.*b22)*rho_fluid*w^2;
beta = (b11 + kFz./rho_fluid/w.^2.*b21).*kFz;
gamma = kFz/rho_fluid/w^2.*b21 + b22;
delta = b22 - kFz/rho_fluid/w^2.*b21;
Rnew = (alpha - beta)./(alpha + beta);
Tnew = delta - gamma.*(alpha - beta)./(alpha + beta);

%% Plot the transmission and reflection
figure
ax(1) = subplot(211);
plot(theta, abs(R), '.', theta, abs(Rnew), 'o')
set(gca, 'YScale', 'log')
ax(2) = subplot(212);
plot(theta, angle(R), '.', theta, angle(Rnew), 'o')

figure
ax(1) = subplot(211);
plot(theta, abs(T), '.', theta, abs(Tnew), 'o')
set(gca, 'YScale', 'log')
ax(2) = subplot(212);
plot(theta, angle(T), '.', theta, angle(Tnew), 'o')

%% Plot the alpha and beta
% figure
% subplot(211)
% plot(theta, abs(alpha), theta, abs(beta))
% subplot(212)
% plot(theta, angle(alpha), theta, angle(beta))

% %% Plot B11
% figure
% ax(3) = gca;
% plot(theta, a22, theta, a21.*a42./a41)
% title('B11')
% figure
% ax(4) = gca;
% plot(theta, b11)
% title('B11')
% 
% %% Plot B22
% figure
% plot(theta, a23, theta, a21.*a43./a41);
% title('B22')
% figure
% plot(theta, b22)
% title('B22')
% 
% %% Plot B12
% figure
% plot(theta, a23, theta, a21.*a43./a41);
% title('B12')
% figure
% plot(theta, b12)
% title('B12')
% 
% %% Plot B21
% figure
% plot(theta, a32, theta, a31.*a42./a41)
% title('B21')
% figure
% plot(theta, b21)
% title('B21')
% 
% %% Determinant
% determinant = arrayfun(@(x11, x12, x21, x22) det([x11, x12; x21, x22]), b11, b12, b21, b22);
% figure
% plot(theta, determinant)
% xlabel('theta')
% ylabel('det(B)')
% 
% %% Compare B12 and B22
% figure
% semilogy(theta, abs(real(b12)))
% hold all
% semilogy(theta, abs(real(b22)))
% legend('abs(Real(B12))', 'abs(Real(B22))')
