% asm/integrands
% Version 100 13-oct-2014
%
% Functions that compute the integrands.
%
% integral* functions are called by the computeAsmIntegral and handles the
% model parameters, and then performes the integral for all the
% frequencies.
%
% Parameter handleing functions:
% integralFluidSolidFluid - Computes the integrand for an axial symmetric
% model using integrandFluidSolidFluidAxialSymmetric. Handles point
% sources, plane pistons and focused transducers. The plate model can be a
% perfect reflector, fluid-solid-fluid model transmission or reflection,
% where the fluids on each side of the solid are identical.
%
%
% Integrands
% integrandFluidSolidFluidAxialSymmetric - 
%
%
%