function [result, p] = parseAsmInput(varargin)
% [result, p] = parseAsmInput() returns the default parameters in
% the struct, result, and a inputParser object, p.
% 
% [result, p] = parseAsmInput('param1', value1, 'param2', value2, ...) returns 
% the parameter struct, result, with the parameters given in the argument, and
% default parameter values for the rest.
%
% 
% Parameters:
%
% Properties of the solid:
% Name - Description - Default value - Unit
% thickness   - Solid plate thickness         - 10e-3 - m
% cs          - Shear speed of sound          - 3158  - m/s
% cp          - Compressional speed of sound  - 5850  - m/s
% rho_solid   - Density of the solid          - 7850  - m/s
% alphaLambda - Absorption times wavelength   - 0     - dB
%                in the solid
% 
% Propertios of the fluid
% Name - Description - Default value - Unit
% cf        - Speed of sound - 342.21 - m/s
% rho_fluid - Density        -   1.5  - kg/m^3
%
% Geometry of the setup
% Name - Description - Default value - Unit
% aTx  - Radius of the transmitter - 9e-3 - m
% aRx  - Radius of the receiver    - 3e-3 - m
% distanceTx - Distance from transmitter to target - 40e-3 - m
% distanceRx - Distance from receiver to target    - 40e-3 - m
% displaceRx - Displacement of the receivers       -  0    - m
%               from the acoustical axis
% tx_focus    - Focal distance of the transmitter   - []    - m
% reflection  - If true, the reflection coefficient - true  - NA
%                is used, otherwise the 
%                transmission coefficient is used
% perfectReflection - If true, the reflection       - false - NA
%                      coefficient is set to 1.
%                      Overrides the reflection
%                      parameter.
% gamma  - Misalignment between transmitter and - 0    - radians
%                the solid plate
%
% Sampling and integration parameters
% f - Frequency to compute the integral - (0:2^15-1)*2e6/2^15 - Hz
% thetamax - Upper integration limit for $\theta_z$ - 0.8 - radians
% thetamin - Lower integration limit for $\theta_z$ - 0   - radians



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
def_perfectReflection = false;
def_gamma = 0;
def_focusTx = 0; % Focal distance (m)
def_focusRx = 0; % Focal distance (m)

% Sampling stuff
def_f = (0:2^15-1)*2e6/2^15;
def_thetamax = 0.8;
def_thetamin = 0;

% Version control
gitInfo = struct('branch', '', 'hash', '', 'remote', '', 'url', '');


%% Parse input
p = inputParser();


% Validator, true if x is a scalar greater than zero (gtz)
validatorsgtz = @(x) validateattributes(x,...
    {'numeric'}, {'scalar', 'real', 'nonempty', '>', 0});

% Validator, true if x is a scalar greater or equal to zero (gtz)
validatorsgoretz = @(x) validateattributes(x,...
    {'numeric'}, {'scalar', 'real', 'nonempty', '>=', 0});

% Validator, true if x is a vector of real numbers
validatorfunvec = @(x) validateattributes(x,...
    {'numeric'}, {'vector', 'real', 'nonempty'});

% Validator, true if x is a scalar and a real number
validatorfuncreal = @(x) validateattributes(x,...
    {'numeric'}, {'scalar', 'real', 'nonempty'});


% Solid
p.addParameter('thickness', def_thickness, validatorsgtz)
p.addParameter('cp', def_cp, validatorsgtz);
p.addParameter('cs', def_cs, validatorsgtz);
p.addParameter('rho_solid', def_rho_solid, validatorsgtz);
p.addParameter('alphaLambda_dB', def_alphaLambda, validatorsgoretz);


% Fluid
p.addParameter('cf', def_cf, validatorsgtz);
p.addParameter('rho_fluid', def_rho_fluid, validatorsgtz);


% Geometric setup
p.addParameter('aTx', def_aTx, validatorsgtz);
p.addParameter('aRx', def_aRx, validatorsgtz);
p.addParameter('distanceTx', def_distanceTx, validatorsgtz);
p.addParameter('distanceRx', def_distanceRx, validatorsgtz);
p.addParameter('displaceRx', def_displaceRx, validatorsgoretz);
p.addParameter('reflection', def_reflection, @islogical);
p.addParameter('perfectReflection', def_perfectReflection, @islogical);
p.addParameter('gamma', def_gamma, validatorfuncreal);
p.addParameter('focusTx', def_focusTx, validatorsgoretz);
p.addParameter('focusRx', def_focusRx, validatorsgoretz);


% Sampling stuff
p.addParameter('f', def_f, validatorfunvec);
p.addParameter('thetamax', def_thetamax, validatorsgtz);
p.addParameter('thetamin', def_thetamin, validatorsgoretz);


% Admin stuff
p.addParameter('filenamevars', {}, @iscell);
p.addParameter('savemat', true, @islogical);
p.addParameter('debug', false, @islogical);


% Version control
p.addParameter('gitInfo', gitInfo, @isstruct);


% Info to screen
p.addParameter('progressbar', false, @islogical);


% Parse and validate filenamevars
p.parse(varargin{:});
validateFilenameVars(p);
result = p.Results;
end


function t = validateFilenameVars(p)
fnv = p.Results.filenamevars;
idf = cellfun(@(x) any(strcmp(x, p.Parameters)), fnv);
if nnz(idf) ~= length(fnv)
    t = false;
    error('HW:INPUTERROR', 'Unknown filenamevar %s. ', fnv{~idf});
end
t = true;
end
