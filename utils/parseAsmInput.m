function [result, p] = parseAsmInput(varargin)
%% Default values
% Solid
def_thickness = 10e-3; % m
def_cs = 3158; % m/s
def_cp = 5850; % m/s
def_rho_solid = 7850; % m/s
def_alphaLambda = 0.008; % dB
% Fluid
def_cf = 342.21; % m/s
def_rho_fluid = 1.5; % kg/m^3
% Geometric setup
def_aTx = 9e-3; % Radius (m)
def_aRx = 3e-3; % Radius (m)
def_distanceTx = 40e-3; % m
def_distanceRx = 40e-3; % m
% Sampling stuff
def_f = (0:2^15-1)*2e6/2^15;
def_thetamax = 0.8;

%% Parse input
p = inputParser();

% Validator, true if x is a scalar greater than zero (gtz)
validatorfun = @(x) validateattributes(x,...
    {'numeric'}, {'scalar', 'real', 'nonempty', '>', 0});
validatorfunvec = @(x) validateattributes(x,...
    {'numeric'}, {'vector', 'real', 'nonempty'});
% Solid
p.addParamValue('thickness', def_thickness, validatorfun)
p.addParamValue('cp', def_cp, validatorfun);
p.addParamValue('cs', def_cs, validatorfun);
p.addParamValue('rho_solid', def_rho_solid, validatorfun);
p.addParamValue('alphaLambda_dB', def_alphaLambda, validatorfun);
% Fluid
p.addParamValue('cf', def_cf, validatorfun);
p.addParamValue('rho_fluid', def_rho_fluid, validatorfun);
% Geometric setup
p.addParamValue('aTx', def_aTx, validatorfun);
p.addParamValue('aRx', def_aRx, validatorfun);
p.addParamValue('distanceTx', def_distanceTx, validatorfun);
p.addParamValue('distanceRx', def_distanceRx, validatorfun);
% Sampling stuff
p.addParamValue('f', def_f, validatorfunvec);
p.addParamValue('thetamax', def_thetamax, validatorfun);
% Admin stuff
p.addParamValue('filenamevars', {}, @iscell);
p.addParamValue('savemat', false, @islogical);
p.addParamValue('debug', true, @islogical);

% Parse
p.parse(varargin{:});
result = p.Results;
end