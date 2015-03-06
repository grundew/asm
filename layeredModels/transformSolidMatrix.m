function B = transformSolidMatrix(A)
% B = transformSolidMatrix(A)
%
%
% Function for calculating the solid layer matrix.
% Equation 24 pp. 1582 [1].
% 
% Input:
% A - 4x4 Solid matrix
% 
%
% Output:
% B - 2x2 Matrix
%
% [1] - A new efficient algorithm to compute the exact reflection
%       and transmission factors for plane waves in 
%       layered absorbing media (liquids and solids)
%       Cervenka, Pierre and Challande, Pascal,
%
%       The Journal of the Acoustical Society of America,
%       89, 1579-1589 (1991), DOI:http://dx.doi.org/10.1121/1.400993



B(1, 1) = A(2, 2) - A(2, 1)*A(4, 2)/A(4, 1);
B(1, 2) = A(2, 3) - A(2, 1)*A(4, 3)/A(4, 1);
B(2, 1) = A(3, 2) - A(3, 1)*A(4, 2)/A(4, 1);
B(2, 2) = A(3, 3) - A(3, 1)*A(4, 3)/A(4, 1);
end