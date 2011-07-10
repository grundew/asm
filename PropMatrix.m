function T = PropMatrix(alpha, beta, z)
% Function propagates the waves with vertical wave number alpha and beta a
% distance z.

T = diag([exp(1i*alpha*z), exp(-1i*alpha*z), exp(1i*beta*z), exp(-1i*beta*z)]);
end