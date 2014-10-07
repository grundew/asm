function p = generateAsmConfig(fluidstr, solidstr)
%p = generateAsmConfig(fluidstr, solidstr) spits out a predefined
% configuration for ASM
% 
% Input:
% str - A string: 'water-alu-water', 'air-steel-air', 'water-steel-air'
%       'water-plexi-water'
% Output
% p - ASM configuration struct

p = parseAsmInput();

switch lower(fluidstr)
    case 'water'
        p.cf = 1500;
        p.rho_fluid = 1000;
    case 'air'
        p.cf = 342;
        p.rho_fluid = 1.5;
end
       
switch lower(solidstr)
    case 'steel'
        p.cp = 5950;
        p.rho_solid = 7950;
        p.cs = 3230;
    case 'alu'
        p.cp = 6400;
        p.cs = 3100;
        p.rho_solid = 2700;
    case 'plexi'
        p.cp = 2400;
        p.cs = 1200;
        p.rho_solid = 1200;
       
    
end

