function [result, p] = parseAsmInput(varargin)
%% Default values

% Solid
def_thickness = 10e-3; % m
def_cs = 3158; % m/s
def_cp = 5850; % m/s
def_rho_solid = 7850; % m/s
def_alphaLambda = 0; % dB

% Fluid
def_cf = 342.21; % m/s
def_rho_fluid = 1.5; % kg/m^3

% Geometric setup
def_aTx = 9e-3; % Radius (m)
def_aRx = 3e-3; % Radius (m)
def_distanceTx = 40e-3; % m
def_distanceRx = 40e-3; % m
def_displaceRx = 0; % m;
def_reflection = true;
def_alpha_plate = 0;
def_tx_focus = false;

% Sampling stuff
def_f = (0:2^15-1)*2e6/2^15;
def_thetamax = 0.8;

%% Parse input
p = inputParser();

% Validator, true if x is a scalar greater than zero (gtz)
validatorsgtz = @(x) validateattributes(x,...
    {'numeric'}, {'scalar', 'real', 'nonempty', '>', 0});
validatorsgoretz = @(x) validateattributes(x,...
    {'numeric'}, {'scalar', 'real', 'nonempty', '>=', 0});
validatorfunvec = @(x) validateattributes(x,...
    {'numeric'}, {'vector', 'real', 'nonempty'});
validatorfuncreal = @(x) validateattributes(x,...
    {'numeric'}, {'scalar', 'real', 'nonempty'});

% Solid
p.addParamValue('thickness', def_thickness, validatorsgtz)
p.addParamValue('cp', def_cp, validatorsgtz);
p.addParamValue('cs', def_cs, validatorsgtz);
p.addParamValue('rho_solid', def_rho_solid, validatorsgtz);
p.addParamValue('alphaLambda_dB', def_alphaLambda, validatorsgoretz);

% Fluid
p.addParamValue('cf', def_cf, validatorsgtz);
p.addParamValue('rho_fluid', def_rho_fluid, validatorsgtz);

% Geometric setup
p.addParamValue('aTx', def_aTx, validatorsgtz);
p.addParamValue('aRx', def_aRx, validatorsgtz);
p.addParamValue('distanceTx', def_distanceTx, validatorsgtz);
p.addParamValue('distanceRx', def_distanceRx, validatorsgtz);
p.addParamValue('displaceRx', def_displaceRx, validatorsgoretz);
p.addParamValue('reflection', def_reflection, @islogical);
p.addParamValue('alpha_plate', def_alpha_plate, validatorfuncreal);
p.addParamValue('tx_focus', def_tx_focus, @islogical);

% Sampling stuff
p.addParamValue('f', def_f, validatorfunvec);
p.addParamValue('thetamax', def_thetamax, validatorsgtz);

% Admin stuff
p.addParamValue('filenamevars', {}, @iscell);
p.addParamValue('savemat', true, @islogical);
p.addParamValue('debug', false, @islogical);

% Parse
p.parse(varargin{:});
result = p.Results;
end