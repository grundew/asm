
nt = 100;
thetamax = 0.25;
theta_F = linspace(0, thetamax, nt);
d = 10e-3;

nf = 3000;
v = 5850;
vShear = 3100;
fres = v/d/2;

% Setup the model
f =  linspace(0.9*fres, 1.05*fres, nf)';
fluid1 = struct('v', 1500, 'density', 1000);
fluid3 = fluid1;
layer = struct('v', v, 'density', 7850, 'vShear', vShear);

% Fluid-solid-fluid model
model = MultiLayerModel(fluid1, layer, fluid3, d);

% L is half the thickness
L = 0.5*model.thickness;
rho_F = model.fluid.density;
rho_S = model.solid.density;
c_L = model.solid.v;
c_S = model.solid.vShear;
c_F = model.fluid.v;

% Calculate angles (independent of frequency)
theta_L = asin(c_L/c_F*sin(theta_F));
theta_S = asin(c_S/c_F*sin(theta_F));

% Calculate wave vectors in the solid, total and vertical
% Total (angle independent)
h = 2*pi*f./c_L;
k = 2*pi*f./c_S;
% Vertical (angle and frequency dependent)
hz = bsxfun(@(x, y) x.*cos(y), h, theta_L);
kz = bsxfun(@(x, y) x.*cos(y), k, theta_S);
% Horizontal (angle and frequency dependent)
K = bsxfun(@(x, y) x.*sin(y), 2*pi*f/c_F, theta_F);

% Normalized horizontal wavelength to thickness
% This is equal to 2*pi at resonance?
alpha_L = hz*L;
alpha_S = kz*L;

%% Helpfull factors
% Frequency independent
DS = rho_S*c_S/(rho_F*c_F)*cos(theta_F)./cos(theta_S);
DL = rho_S*c_L/(rho_F*c_F)*cos(theta_F)./cos(theta_L);
A = c_L/c_S*cos(theta_S)/cos(theta_L);
% Frequency dependent
EL = bsxfun(@(x, y) cos(2*y).^2./sin(2*x), alpha_L, theta_S);
ES = bsxfun(@(x, y) sin(2*y).^2./sin(2*x), alpha_S, theta_S);

%% Calculate reflection and transmission coefficients via N and M
N1 = bsxfun(@(x, y) x.*y, DL, EL);
N2 = bsxfun(@(x, y) x.*y, DS, ES);
M1 = bsxfun(@(x, y) x.*y, DL, EL.*cos(2*alpha_L));
M2 = bsxfun(@(x, y) x.*y, DS, ES.*cos(2*alpha_S));
N = N1 + N2;
M = M1 + M2;

% Transmission and reflection
T = 2*N./(2*M + 1i*(M.^2 - N.^2 - 1));
R = 1i*(M.^2 - N.^2 + 1)./(2*M + 1i*(M.^2 - N.^2 - 1));
T(f==0, :) = ones(nnz(f==0), length(theta_F));
R(f==0, :) = zeros(nnz(f==0), length(theta_F));

%% Plot the parameters
x = f/fres;
plotit = @(z) imagesc(x, theta_F, abs(z).');

% Plot reflection coefficient
figure
plotit(R), axis xy, colormap gray

% Plot the sum and difference of alpha_X square
% Can drop the sum?
figure
alphasum = (alpha_L + alpha_S);
alphadiff = (alpha_L - alpha_S);
subplot(311),plotit(1 - cos(alphasum).^2),colorbar,axis xy,title('1 - cos(alpha_L + alpha_S).^2')
subplot(312),plotit(1 - cos(alphadiff).^2),colorbar,axis xy,title('1 - cos(alpha_L - alpha_S).^2')
subplot(313),plotit(1 - cos(alphadiff).^2 - cos(alphadiff).^2),colorbar,axis xy,title('1-alphasum-alphadiff')

%% Calculate the three terms in the numerator of R
% Scalars
AL = rho_S*c_L/rho_F/c_F;
AS = rho_S*c_S/rho_F/c_F;
% Theta dependent
BL = cos(theta_F)./cos(theta_L);
BS = cos(theta_F)./cos(theta_S);
% Frequency dependent
DL = bsxfun(@rdivide, cos(2*theta_S).^2, sin(2*alpha_L));
DS = bsxfun(@rdivide, sin(2*theta_S).^2, sin(2*alpha_S));

% Gather terms. M^2 - N^2 + 1 = -term1 - term2 - term3 +1; 
term1 = AL^2*bsxfun(@times, DL.*sin(2*alpha_L), BL).^2;
term2 = AS^2*bsxfun(@times, DS.*sin(2*alpha_S), BS).^2;
term3 = 2*AL*AS*bsxfun(@times,...
    DL.*DS.*(1 - cos(2*alpha_L).*cos(2*alpha_S)),...
    BL.*BS);

% Check the result. Seems correct now to within 1e-9
Rcheck = 1i*(-term1 - term2 - term3 + 1)./(2*M + 1i*(M.^2 - N.^2 - 1));
% Plot the new reflection coefficient to check
figure
plotit(Rcheck), axis xy, colormap gray,title('Reflection coefficient check')

% Plot the three terms
figure
plotit(term1), axis xy, title term1
figure
plotit(term2), axis xy, title term2
figure
plotit(term3), axis xy, title term3


