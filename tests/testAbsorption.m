function testAbsorption()

nfft = 2^12;
fs = 5e6;
t = (0:nfft-1)/fs;
f = (0:nfft-1)*fs/nfft;
c = 1500;
f0 = 200e3;
alphaLambda_dB = 0.008;

% Pulse
t_ex = (0:1/fs:50e-6);
x_ex = conj(hilbert(sin(2*pi*f0*t_ex)));
X_ex = ifft(x_ex, nfft);

% Propagate wave one wavelength
lambda = c/f0;

% With damping (complex k_z)
alpha_L = alphaLambda_dB*log(10)/10*f0/c;
c_L = c./(1 + 1i*alpha_L*c./f/2/pi);
c_L(f==0) = c;
k_z = 2*pi*f./c_L;


%% Propagate the wave
figure
plot(t_ex, real(x_ex))
hold on
prop = exp(1i*k_z*lambda);
x_prop = fft(X_ex.*prop, nfft);

plot(t, real(x_prop))

a_0 = max(real(x_ex));
a_lambda = max(real(x_prop));
err = abs(10*log10(a_0/a_lambda)-alphaLambda_dB);
assert(err < 10^-6,...
    'Calculated loss is not equal to the input')
end