function kappa = poisson2kappa(mu)
% mu = kappa2poission(kappa)
% 
% Calculates Poisson's ratio, $\mu$ from $\kappa = c_p/c_s$.
%
% Input:
% 
% mu - $\mu$ Poisson's ratio (scalar, vector or matrix)
% 
% Output:
%
% kappa - $\kappa$ = c_p/c_s$


kappa = sqrt(2*(1 - mu)./(1 - 2*mu));


end