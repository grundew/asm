function mu = kappa2poisson(kappa)
% mu = kappa2poission(kappa)
% 
% Calculates Poisson's ratio, $\mu$ from $\kappa = c_p/c_s$.
%
% Input:
% 
% kappa - $\kappa$ = c_p/c_s$ (scalar, vector or matrix)
% 
% Output:
% mu - $\mu$ Poisson's ratio


mu = 0.5*(kappa.^2 - 2)./(kappa.^2 - 1);


end