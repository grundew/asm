clear
fluid = materials.Fluid;
fluid.v = 1500;
fluid.density = 1000;
solid = materials.LinearElastic;
solid.v = 5900;
solid.vShear = 3150;
solid.density = 7850;
d = 12.4e-3;
model = MultiLayerModel(fluid, solid, fluid, d);

%f = 0:0.5e3:600e3;
%lambda = solid.v./f;

theta_c = asin(fluid.v/solid.vShear);
%theta_c = asin(fluid.v/solid.v);

theta = 0.01:0.01:theta_c;
%theta = 0;
%theta = 0.08;

% First shear
xl = 0.52;
xh = 0.64;

% Second compressional
%xl = 1.8;
%xh = 2.4;

%nn = 2^10;
nn = 2.^13;

figure
for n = nn
    
    x = linspace(0, 3, n);
    f = 0.5*x*solid.v/d;
    [V, W, k_hor, k_vert_L, debug] = fluidSolidFluid(f, theta, model);
    %     figure
    %     ax(1) = subplot(211);
    %     plot(x, real(V), '.');
    %     ax(2) = subplot(212);
    %     plot(x, imag(V), '.');
    %     linkaxes(ax, 'x');
    
    % Estimate energy in shear resonance
    id = find(x > xl & x < xh);
    energy = sum(abs(V(id, :)).^2)/length(id);
    hold all
    plot(theta, energy)

end

legend(cellfun(@(x) num2str(x), num2cell(nn), 'UniformOutput', 0))
xlabel('\theta')
ylabel('Energy')
figure
subplot(211)
plot(x, real(V), '.')
subplot(212)
plot(x, imag(V), '.')