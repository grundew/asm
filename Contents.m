% asm
% Version 100 13-oct-2014
%
%
% computeAsmIntegral() is the main function that integrates the specified
% function wrt the $\theta_z$ angle. It takes a function handle and a
% parameter struct as input.
%
% Additional parameter-value inputs are sent to the integration function,
% quadgk. See the online MATLAB documentation for quadgk for details.
% 
% The functions in the <matlab:help('intergands/Content') integrals> 
% subfolder defines the function that is integrated in computeAsmIntegral.
%
%
%
% Dependencies:
% Version control
% The following third party apps must be included in the MATLAB path as
% well as git.exe in the windows path. This is optional, but when available
% the git revision information will be stored with the computed data.
%
% <https://github.com/manur/MATLAB-git MATLAB-git>
% <https://github.com/quantentunnels/matlab-getgitinfo matlab-getgitinfo>
%
%
% Progress bar (ASCII):
% <http://www.mathworks.com/matlabcentral/fileexchange/8564-progress progress>
%
%
%
% Example:
% p = generateAsmConfig('air', 'steel', 'distanceRx', 10e-2);
% func = @integralFluidSolidFluid;
% [V, f] = computeAsmIntegral(func, p, 'MaxIntervalCount', 2000);
%
%
%
% See also quadgk, generateAsmConfig, computeAsmIntegral,
% asm/integrals
%
%
%
% Unit tests:
% Run all unit tests with the command:
%   runasmtests();
%
%
%
% References:
% 1. Orofino 1993 - http://dx.doi.org/10.1121/1.405408
%
%
%
% Files:
%   computeAsmIntegral  - Main function. Computes the integral over $\theta_z$.
%   planeWaveTimeSignal - Computes time signal for a single plane wave
%                         model
%   runasmtests         - Runs all test cases in the tests folder
%
%
% Folders:
% 	integrands			- Angular spectrum integrands, passed as function handles
%						  to the computeAsmIntegral function.
%	layeredModels		- Functions for calculating reflection/transmission
%						  coefficients for plane, layered media.
%	semiInfiniteMedia	- Functions for calculating reflection/transmission 
%						  coefficients for the interface between
%						  two semi-infinite media
%	utils				- Various functions for handling parameters, batch simulation etc.
%	tests				- Unit tests. Run all tests with runasmtests-command.
%	deprecated			- Outdated functions
%
%
%
