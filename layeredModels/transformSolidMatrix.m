function B = transformSolidMatrix(A)
% Function for calculating the solid layer matrix.
% 
% B = transformSolidMatrix(A)
%
% J. AcoustS. oc.A m.8 9 (4), Pt. 1, April 1991
% A new efficient algorithm to compute the exact reflection
% and transmission factors for plane waves in layered absorbing
% media (liquids and solids)
% Pierre Cervenka and Pascal Challande
% Equation 24 pp. 1582

B(1, 1) = A(2, 2) - A(2, 1)*A(4, 2)/A(4, 1);
B(1, 2) = A(2, 3) - A(2, 1)*A(4, 3)/A(4, 1);
B(2, 1) = A(3, 2) - A(3, 1)*A(4, 2)/A(4, 1);
B(2, 2) = A(3, 3) - A(3, 1)*A(4, 3)/A(4, 1);
end