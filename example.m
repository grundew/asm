fluid = materials.Fluid;
fluid.v = 1450;
fluid.density = 1000;
solid = materials.LinearElastic;
solid.v = 5900;
solid.vShear = 3150;
solid.density = 7850;
d = 12.4e-3;
model = MultiLayerModel(fluid, solid, fluid, d);

f = 100e3:1e3:600e3;
theta = 0.01:0.01:0.1;
[V, W, k_hor, k_vert_L, debug] = fluidSolidFluid(f, theta, model);
figure(2)
subplot(211)
plot(2*pi*f, real(V));
subplot(212)
plot(2*pi*f, imag(V));