function L = Lmatrix(alpha, beta, gamma, mu)
% Function computes the L matrix according to equation 4.3.5 p 99
% Brekhovskihk

L(1, :) = [1i*xi, 1i*xi, -1i*beta, 1i*beta];
L(2, :) = [1i*alpha, -1i*alpha, 1i*xi, 1i*xi];
L(3, :) = [2*mu*xi*gamma, 2*mu*xi*gamma, -2*mu*xi*beta, 2*mu*xi*beta];
L(4, :) = [-2*mu*xi*alpha, 2*mu*xi*alpha, -2*mu*xi*gamma, -2*mu*xi*gamma];
end