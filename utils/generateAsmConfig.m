function result = generateAsmConfig(fluidstr, solidstr, varargin)
%p = generateAsmConfig(fluidstr, solidstr) spits out a predefined
% configuration for ASM
% 
% Input:
% fluidstr - A string: 'water' or 'air'.
% solidstr - A string: 'steel', 'alu' or 'plexi'
% varargin - Additional parameters to parseAsmInput
%
% Output:
% result - ASM configuration struct

[result, p] = parseAsmInput(varargin{:});

s = {'cf', 'cp', 'cs', 'rho_solid', 'rho_fluid'};
if ~all(cellfun(@(x) any(strcmp(p.UsingDefaults, x)), s))
    warning('Some of the parameters specified will be changed to some other default')
end

switch lower(fluidstr)
    
    case 'water'
        result.cf = 1500;
        result.rho_fluid = 1000;
    case 'air'
        result.cf = 342;
        result.rho_fluid = 1.5;
end
       
switch lower(solidstr)
    case 'steel'
        result.cp = 5950;
        result.rho_solid = 7950;
        result.cs = 3230;
    case 'alu'
        result.cp = 6400;
        result.cs = 3100;
        result.rho_solid = 2700;
    case 'plexi'
        result.cp = 2400;
        result.cs = 1200;
        result.rho_solid = 1200;
       
    
end