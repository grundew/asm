%% Frequency, speed of sound and wavelengths
f = 1000e3;
w = 2*pi*f;
cL = 5850;
lambdaL = cL/f;
theta = 0.01;

%% Define damping coefficient
alphaL = 0.008/lambdaL/8.69;

%% Define complex velocity
rL = alphaL*lambdaL/2/pi;
c = cL/(1 + 1i*rL);

%% Compute wave-numbers and check the imaginary part vs alpha
k = w/c;
kz = k*cos(theta);
kr = k*sin(theta);

%% Compare
fprintf('theta < thetaC\n')
fprintf('alphaL: %f\n', alphaL)
fprintf('Im(k): %f\n', imag(k))
fprintf('Im(kz): %f\n', imag(kz))
fprintf('Im(kr): %f\n', imag(kr))

%% Compute wave-numbers at angles above the critical angle
cF = 1500;
thetac = asin(cF/cL);
theta = thetac*1.2;
thetaL = asin(c/cF*sin(theta));

k = w/c;
kzc = k*cos(thetaL);
krc = k*sin(thetaL);

% Compare
fprintf('theta > thetaC\n')
fprintf('alphaL: %f\n', alphaL)
fprintf('Im(k):  %f\n', imag(k))
fprintf('Im(kz): %f\n', imag(kz))
fprintf('Im(kr): %f\n', imag(kr))

% Compute wave numbers without the complex wave speed
thetaL = asin(cL/cF*sin(theta));
k = w/cL;
kz = k*cos(thetaL);
kr = k*sin(thetaL);

fprintf('Whithout the complex wave speed\n')
fprintf('alphaL: %f\n', alphaL)
fprintf('Im(k):  %f\n', imag(k))
fprintf('Im(kz): %f\n', imag(kz))
fprintf('Im(kr): %f\n', imag(kr))