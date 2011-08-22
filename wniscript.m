% Step 1:
% Calculate input matrix
input = inputMatrix(rho, w, k_z);

% Step 2:
% Calculate the matrices related to each layer
a = layerMatrix(rho, w, k_z_S, k_z_L, K, d);

% Step 3:
% Calculate the output matrix
output = outputMatrix(rho, w, k_z_L);

% Step 4:
% Calculate the products of all the similar contiguous matrices.

% A = a(end, :, :);
% for i = length(a)-1:-1:1
%     A = a(i+1, :, :)*a(i, :, :);
% end

% Single layer
A = a;

% Step 5:
% Reduceing of the liquid-solid-liquid sandwiches.
B = transformSolidMatrix(A);