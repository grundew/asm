% asm/integrands
% Version 100 13-oct-2014
%
% Functions that compute the integrands.
%
% integral* functions are called by the computeAsmIntegral and handles the
%           model parameters, and then performes the integral for all the
%           frequencies. The integrands are defined in the integrand*
%           functions.
%
%
% Parameter handling functions:
% integralFluidSolidFluid - Computes the integrand for an axial symmetric
%                           model using
%                           integrandFluidSolidFluidAxialSymmetric. Handles
%                           point sources, plane pistons and focused
%                           transducers. The plate model can be a perfect
%                           reflector, fluid-solid-fluid model transmission
%                           or reflection, where the fluids on each side of
%                           the solid are identical.
%
%
% Integrands:
% integrandFluidSolidFluidAxialSymmetric:
%             Integrand for the axial symmetric model, used by
%             integralFluidSolidFluid. Plate is either embedded in a fluid
%             or a perfect reflector.
%
% integrandFluidSolidFluid2Dangle
%             Integrand for 2D model with misalignment between plate and
%             transducers. Only plane piston. Note that this function
%             should be integrated from $-pi/2$ to $pi/2$.
%
% integrandFluidSolidFluid_withAngle
%             3D model with misaligned plate. Not implemented properly.
%
%
% Sources:
% focusedSourceASM - Angular spectrum of a focused bowl transducer.