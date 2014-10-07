function X = integrandFluidSolidFluid_planepiston(params, varargin)
% I = orofinoIntegrand(theta, f, aRx, aTx,...
%                      c, rho, d1, d3, model,...
%                      alphaLambda_dB, reflection, alpha_plate)
%   orofinoIntegrand is the integrand in equation (28) in Ref. 1. This is
%   used as an input to the waveNumberIntegration function. It uses an
%   analytical expression for the transmission coefficient, assuming that
%   the fluid on both sides of the solid plate have the same properties.
% 
% Input:
% theta_z - Angle (scalar or vector)
% f - Frequency (scalar)
% aRx - Radius of the receiver
% aTx - Raidus of the transmitter
% c - Speed of sound in the propagation fluid
% rho - Density of the propagation fluid
% d1 - Distance from transmitter to the solid plate
% d3 - Distance from receiver to the solid plate
% x0 - Displacement of receiver from the acoustical axis of the transmitter
% model - MultiLayerModel object
% alphaLambda - Damping factor
% reflection - (Boolean) True if reflection coeffecient should be used
% alpha_plate - Angle of misalignment of the plate
%
% Output:
% I - Integrand evaluated at theta
%
% References:
% 1. Orofino, 1992. http://dx.doi.org/10.1121/1.405408

%% Check sanity of parameters and give warnings or errors


%% Samplings stuff
f = params.f;
thetamax = params.thetamax;
thetamin = 0;

%% Unpack parameters
aRx = params.aRx;
aTx = params.aTx;
v_fluid = params.cf;
rho_fluid = params.rho_fluid;
cp = params.cp;
d1 = params.distanceTx;
d3 = params.distanceRx;
alphaLambda_dB = params.alphaLambda_dB;
fres = 0.5*params.cp/params.thickness; %#ok<*NASGU>
x0 = params.displaceRx;
refl = params.reflection;


%% Pack the parameters for the model
fluid1 = struct('v', params.cf, 'density', params.rho_fluid);
fluid3 = fluid1;
layer = struct('v', params.cp, 'density', params.rho_solid, 'vShear', params.cs);

% Fluid-solid-fluid model
model = struct('fluid', [fluid1, fluid3],...
    'solid', layer, 'thickness', params.thickness,...
    'alphaLambda', alphaLambda_dB);

X = zeros(size(f));

for i = 1:length(f)
    % Time it
    if i == 1
        fprintf('Started: %s\n', datestr(now, 'dd-mm-yyyy_HHMMSS'));
        tstart = tic();
    end
    
    w = 2*pi*f(i);
    k = w/v_fluid;
    
    % fun = @(theta_z) integrand(...
    %    theta_z, freq, aRx, aTx,...
    %    v_fluid, d1, d3, model, x0, alphaLambda_dB, refl);
    fun = @(theta_z) integrand(theta_z, f(i), k, d1, d3, aTx, aRx, v_fluid,...
        alphaLambda_dB, refl, model, x0);
    X(i) = quadgk(fun, thetamin, thetamax, 'MaxIntervalCount', 5000);
    
    % Time it
    if i == 500
        tme = toc(tstart)*length(f)/i;
        nn = datenum(0, 0, 0, 0, 0, tme);
        nnvec = datevec(nn);
        fprintf('Estimated time of arrival: %.0f HOURS %.0f MIN %.0f SEC\n',...
            nnvec(end-2), nnvec(end-1), ceil(nnvec(end)));
    end

end

end

function I = integrand(theta_z, f, k, d1, d3, aTx, aRx, c,...
    alphaLambda_dB, reflection, model, x0)
%% Compute wave numbers
q = sin(theta_z);
p = sqrt(1-q.^2);

kr = k*q;
kz = k*p;

%% Transmitter spatial sensitivity
% No angle adjustment
xx = kr*aTx;
W = besselj(1, xx)./xx;
W(xx==0) = 0.5;
PhiTx = 2*pi*aTx^2*W;


%% Receiver spatial sensitivity
xx = kr*aRx;
Wrx = besselj(1, xx)./xx;
Wrx(xx==0) = 0.5;
PhiRx = 2*pi*aRx^2*Wrx;


%% Plate response, angular
% Multiply with wave length and convert from dB to linear
% Loss parameter
if alphaLambda_dB > 0
    % log(10)/10 = 0.2303
    alphaL = alphaLambda_dB*0.2303*f/c;
else
    alphaL = 0;
end


%% Reflection/Transmission coefficient
if reflection
    % TODO: Add loss
    Plate = analyticRTFast(f, theta_z, model, alphaL);
else
    Plate = transmissionCoefficientAnalytical(f, q, model, alphaL);
end


%% Displacement factor
if x0 > 0
    dispRx = besselj(0, x0*kr);
else
    dispRx = 1;
end


%% Phase shift from transmitter to plate and from plate to receiver
Phase = exp(1i*kz*(d1 + d3));


%% Assemble integrand
I = Plate.*k.*q.*dispRx.*PhiRx.*PhiTx.*Phase.*k.*p.^2;


end