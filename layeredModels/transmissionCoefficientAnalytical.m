function T = transmissionCoefficientAnalytical(freq, q, model, alpha_L, alpha_S)
% T = transmissionCoefficientAnalytical(freq, q, model, alpha_L, alpha_S)
%
%
% transmissionCoefficientAnalytical computes the transmission coefficient
% for a solid elastic plate embedded in a fluid (same fluid on both sides).
% Equation (1) in [1] is used.
%
%
% Input:
% freq     - Frequency (scalar, Hz)
% q        - Sine of theta, where theta is the angle of incidence in the fluid.
% model    - A MultiLayerModel object.
% alpha_L  - Damping coefficient (optional, scalar, 1/m). Default = 0.
% alpha_S  - Damping coefficient (optinoal, scalar, 1/m). If not given, and
%            alpha_L is given, alpha_S is calculated based on alpha_L. See
%            Equation 11 in [2].
%
% 
% Output:
% T        - Complex transmission coefficient, length(freq) x length(q)
%
%
% [1] - MEASUREMENTS AND 3D SIMULATIONS OF ULTRASONIC DIRECTIVE BEAM
%       TRANSMISSION THROUGH A WATER-IMMERSED STEEL PLATE
%       Kjetil Daae Lohne, Per Lunde, Magne Vestrheim
%       Proceedings of the 34th Scandinavian Symposium on Physical Acoustics, Geilo 30 January ? 2 February, 2011
%
% [2] - Modal resonance analysis of acoustic transmission and reflection
%       losses in viscoelastic plates
%       Walter Madigosky and Ralph Fiorito
%       JASA 65 page 1105

if length(freq) == 1 && freq == 0
    T = ones(size(q));
    return;
end

% L is half the thickness
L = 0.5*model.thickness;
rho_F = model.fluid.density;
rho_S = model.solid.density;

% The real wave velocities
c_Lr = model.solid.v;
c_Sr = model.solid.vShear;
c_F = model.fluid.v;

%% Handle absorption as complex wave velocity
% Check for complex speeds of sound
if nargin < 4
    alpha_L = 0;
    alpha_S = 0;
elseif nargin < 5
    % Ref(2) equation 11
    alpha_S = (alpha_L/2)*(c_Lr/c_Sr)^3;
end

if alpha_L == 0
    c_L = c_Lr;
    c_S = c_Sr;
else
    % Equation 5 & 6 from Ref(2)
    c_L = c_Lr./(1 + 1i*alpha_L*c_Lr./freq/2/pi);
    c_S = c_Sr./(1 + 1i*alpha_S*c_Sr./freq/2/pi);
end

% Calculate angles (independent of frequency)
theta_L = asin(c_L/c_F*q);
theta_S = asin(c_S/c_F*q);

% Calculate wave vectors in the solid, total and vertical
% Total (angle independent)
h = 2*pi*freq./c_L;
k = 2*pi*freq./c_S;

% Vertical (angle and frequency dependent)
hz = bsxfun(@(x, y) x.*cos(y), h, theta_L);
kz = bsxfun(@(x, y) x.*cos(y), k, theta_S);

% Normalized horizontal wavelength to thickness
% This is equal to 2*pi at resonance?
alpha_L = hz*L;
alpha_S = kz*L;

%% Helpfull factors
% Frequency independent
qcos = sqrt(1-q.^2);
DS = rho_S*c_S/(rho_F*c_F)*qcos./cos(theta_S);
DL = rho_S*c_L/(rho_F*c_F)*qcos./cos(theta_L);

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
idfreq0 = freq==0;
T(:, idfreq0) = 1;
end