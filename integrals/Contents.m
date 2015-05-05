% asm/integrals
% Version 100 13-oct-2014
%
% Functions that compute the integrals.
%
% integral* functions are passed to the computeAsmIntegral as function handles
%           and handles the model parameters, and then performes the integral
%           for all the frequencies. The integrands are defined in the integrands*
%           functions or in the subfunctions.
%
%
% Parameter handling functions:
% integralFluidSolidFluid - Computes the integral for an axial symmetric
%                           model using
%                           integrandFluidSolidFluidAxialSymmetric. Handles
%                           point sources, plane pistons and focused
%                           transducers. The plate model can be a perfect
%                           reflector, fluid-solid-fluid model transmission
%                           or reflection, where the fluids on each side of
%                           the solid are identical.
% 
% integralFluidSolidFluidMisalignment2D - Computes the integral for an
%                           axial symmetric model with a misalignment of
%                           the solid plate. Handles plane piston
%                           transmitter end receiver. The plate model can
%                           be a perfect reflector or a solid plate
%                           embedded in a fluid.
%
%
% Integrands:
% integrandFluidSolidFluidAxialSymmetric:
%             Integrand for the axial symmetric model, used by
%             integralFluidSolidFluid. Plate is either embedded in a fluid
%             or a perfect reflector.
%
% integrandFluidSolidFluid3Dangle
%             3D model with misaligned plate. Not tested.
%
%
% Sources:
% focusedSourceASM - Angular spectrum of a focused bowl transducer.
