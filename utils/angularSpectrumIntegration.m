function I = angularSpectrumIntegration(f, a, b, h, errTol)
% Use the composite trapezoidal rule to integrate. The function includes
% more and points until the value of the integral stops changing more than
% errTol.

if nargin < 5
    errTol = 1e-6;
end
% Random init number
abserr = errTol + 10;

% Add the average of the end points
I = 0.5*h*(f(a) + f(b));
% Random init number
Iprev = 1e6;

while abserr > errTol
    x = a + (1:n).*h;
    y = f(x);
    I = I + h*sum(y(2:end-1));
    abserr = abs(I - Iprev);
    Iprev = I;
end

end